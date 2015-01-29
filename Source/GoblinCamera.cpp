#include "GoblinCamera.h"
#include "GoblinFilm.h"
#include "GoblinRay.h"
#include "GoblinSampler.h"
#include "GoblinUtils.h"

namespace Goblin{

    Camera::Camera(const Vector3& position, const Quaternion& orientation,
        float fov, float zn, float zf, 
        float lensRadius, float focalDistance, Film* film):
        mPosition(position), mOrientation(orientation),
        mFOV(fov), mZNear(zn), mZFar(zf), 
        mLensRadius(lensRadius), mFocalDistance(focalDistance), 
        mFilm(film) {
        float xRes = static_cast<float>(film->getXResolution());
        float yRes = static_cast<float>(film->getYResolution());
        mAspectRatio = xRes / yRes; 
        mProj = matrixPerspectiveLHD3D(mFOV, mAspectRatio, mZNear, mZFar);
        update();
    }

    //TODO remove this
    Camera::Camera() {
        mPosition = Vector3(0.0f, 0.0f, 0.0f);
        mOrientation = Quaternion::Identity;

        mProj = matrixPerspectiveLHD3D(PI * 0.25f, 1.3333f, 0.5f, 1000.0f);
        mFilm = NULL;
        update();
    }

    Camera::~Camera() {
        if(mFilm != NULL) {
            delete mFilm;
            mFilm = NULL;
        }
    }

    float Camera::generateRay(const Sample& sample, Ray* ray) {
        float xNDC = +2.0f * sample.imageX / mFilm->getXResolution() - 1.0f;
        float yNDC = -2.0f * sample.imageY / mFilm->getYResolution() + 1.0f;
        // from NDC space to view space
        // xView = xNDC * zView * tan(fov / 2) * aspectRatio
        // yView = yNDC * zView * tan(fov / 2)
        // in projection matrix pro,
        // pro[0][0] = 1 / (tan(fov / 2) * aspectRatio)
        // pro[1][1] = 1 / (tan(fov / 2))
        float zView = 1.0f;
        float xView = xNDC / mProj[0][0];
        float yView = yNDC / mProj[1][1];
        Vector3 viewDir = normalize(Vector3(xView, yView, zView));
        if(mLensRadius == 0.0f) {
            // view space to world space
            ray->o = mPosition;
            ray->d = mOrientation * viewDir;
        } else {
            // perturb the ray direction for DOF effect
            float ft = mFocalDistance / viewDir.z;
            Vector3 pFocus = viewDir * ft;
            Vector2 lensSample = mLensRadius *
                uniformSampleDisk(sample.lensU1, sample.lensU2);
            Vector3 viewOrigin(lensSample.x, lensSample.y, 0.0f);
            ray->o = mOrientation * viewOrigin + mPosition;
            ray->d = mOrientation * normalize(pFocus - viewOrigin);
        }
        ray->mint = 0.0f;
        ray->maxt = INFINITY;
        ray->depth = 0;
        return 1.0f;
    }

    const Matrix4& Camera::getWorldMatrix() {
        if(!isUpdated()) {
            update();
        }
        return mWorld;
    }

    const Matrix4& Camera::getViewMatrix() {
        if(!isUpdated()) {
            update();
        }
        return mView;
    }

    const Matrix4& Camera::getProjectionMatrix() {
        return mProj;
    }

    const Vector3 Camera::getLook() const {
        return mOrientation * Vector3::UnitZ;
    }

    const Vector3 Camera::getUp() const {
        return mOrientation * Vector3::UnitY;
    }

    const Vector3 Camera::getRight() const {
        return mOrientation * Vector3::UnitX;
    }

    void Camera::roll(float angle) {
        rotate(mOrientation * Vector3::UnitZ, angle);
    }

    void Camera::pitch(float angle) {
        rotate(mOrientation * Vector3::UnitX, angle);
    }

    void Camera::yaw(float angle) {
        rotate(mOrientation * Vector3::UnitY, angle);
    }

    void Camera::rotateX(float angle) {
        rotate(Vector3::UnitX, angle);
    }

    void Camera::rotateY(float angle) {
        rotate(Vector3::UnitY, angle);
    }

    void Camera::rotateZ(float angle) {
        rotate(Vector3::UnitZ, angle);
    }
    
    void Camera::rotate(const Vector3& axis,float angle) {
        mOrientation = normalize(Quaternion(axis, angle) * mOrientation);
        mIsUpdated = false;
    } 

    void Camera::translate(const Vector3& d) {
        mPosition += d;
        mIsUpdated = false;
    }

    Vector3 Camera::worldToScreen(const Vector3& pWorld) const {
        Vector4 pView = mView * Vector4(pWorld, 1.0f); 
        Vector4 pNDC = mProj * pView;
        pNDC /= pNDC.w;
        float screenX = (pNDC.x + 1.0f) * 0.5f * mFilm->getXResolution();
        float screenY = (1.0f - pNDC.y) * 0.5f * mFilm->getYResolution();
        return Vector3(screenX, screenY, pView.z);
    }
    
    void Camera::update() {
        Matrix4 R = mOrientation.toMatrix();

        Vector3 right(R[0][0], R[1][0], R[2][0]);
        Vector3    up(R[0][1], R[1][1], R[2][1]);
        Vector3  look(R[0][2], R[1][2], R[2][2]);

        mWorld = Matrix4(
             right.x, up.x, look.x, mPosition.x,
             right.y, up.y, look.y, mPosition.y,
             right.z, up.z, look.z, mPosition.z,
                0.0f, 0.0f,   0.0f,        1.0f);

        float x = -dot(mPosition, right);
        float y = -dot(mPosition, up);
        float z = -dot(mPosition, look);

        mView = Matrix4(
            right.x, right.y,  right.z,    x,
               up.x,    up.y,     up.z,    y,
             look.x,  look.y,   look.z,    z,
               0.0f,    0.0f,     0.0f, 1.0f);

        mIsUpdated = true;
    }

    Camera* PerspectiveCameraCreator::create(const ParamSet& params, 
        Film* film) const {
        Vector3 position = params.getVector3("position");
        Quaternion orientation;
        if(params.hasVector3("euler")) {
            Vector3 xyz = params.getVector3("euler", Vector3::Zero);
            orientation = Quaternion(Vector3::UnitZ, radians(xyz.z)) * 
                Quaternion(Vector3::UnitY, radians(xyz.y)) *
                Quaternion(Vector3::UnitX, radians(xyz.x));
        } else {
            Vector4 q = params.getVector4("orientation", Vector4(1, 0, 0, 0));
            orientation = Quaternion(q[0], q[1], q[2], q[3]);
        }
        float fov = params.getFloat("fov", 60.0f);
        float zn = params.getFloat("near_plane", 0.1f);
        float zf = params.getFloat("far_plane", 1000.0f);
        float lensRadius = params.getFloat("lens_radius", 0.0f);
        float focalDistance = params.getFloat("focal_distance", 1e+10f);
        return new Camera(position, orientation, radians(fov), 
            zn, zf, lensRadius, focalDistance, film);
    }
}
