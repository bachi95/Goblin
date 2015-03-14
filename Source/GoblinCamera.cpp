#include "GoblinCamera.h"
#include "GoblinFilm.h"
#include "GoblinRay.h"
#include "GoblinSampler.h"
#include "GoblinUtils.h"

namespace Goblin{

    Camera::Camera(const Vector3& position, const Quaternion& orientation,
        float zn, float zf, Film* film):
        mPosition(position), mOrientation(orientation),
        mZNear(zn), mZFar(zf), mFilm(film), mIsUpdated(false) {
        float xRes = static_cast<float>(film->getXResolution());
        float yRes = static_cast<float>(film->getYResolution());
        mAspectRatio = xRes / yRes; 
    }

    Camera::~Camera() {
        if(mFilm != NULL) {
            delete mFilm;
            mFilm = NULL;
        }
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

    PerspectiveCamera::PerspectiveCamera(const Vector3& position, 
        const Quaternion& orientation, float fov, float zn, float zf, 
        float lensRadius, float focalDistance, Film* film):
        Camera(position, orientation, zn, zf, film),
        mFOV(fov), mLensRadius(lensRadius), mFocalDistance(focalDistance) {
        mProj = matrixPerspectiveLHD3D(mFOV, mAspectRatio, mZNear, mZFar);
        update();
    }

    float PerspectiveCamera::generateRay(const Sample& sample, 
        RayDifferential* ray) const {
        float invXRes = mFilm->getInvXResolution();
        float invYRes = mFilm->getInvYResolution();
        float xNDC = +2.0f * sample.imageX * invXRes - 1.0f;
        float yNDC = -2.0f * sample.imageY * invYRes + 1.0f;
        float dxNDC = +2.0f * (sample.imageX + 1.0f) * invXRes - 1.0f;
        float dyNDC = -2.0f * (sample.imageY + 1.0f) * invYRes + 1.0f;
        // from NDC space to view space
        // xView = xNDC * zView * tan(fov / 2) * aspectRatio
        // yView = yNDC * zView * tan(fov / 2)
        // in projection matrix pro,
        // pro[0][0] = 1 / (tan(fov / 2) * aspectRatio)
        // pro[1][1] = 1 / (tan(fov / 2))
        float zView = 1.0f;
        float xView = xNDC / mProj[0][0];
        float yView = yNDC / mProj[1][1];
        Vector3 viewDir = Vector3(xView, yView, zView);

        float dxView = dxNDC / mProj[0][0];
        Vector3 dxViewDir(dxView, yView, zView);
        float dyView = dyNDC / mProj[1][1];
        Vector3 dyViewDir(xView, dyView, zView);

        if(mLensRadius == 0.0f) {
            // view space to world space
            ray->o = ray->dxOrigin = ray->dyOrigin = mPosition;
            ray->d = mOrientation * normalize(viewDir);
            ray->dxDir = mOrientation * normalize(dxViewDir);
            ray->dyDir = mOrientation * normalize(dyViewDir);
            
        } else {
            // perturb the ray direction for DOF effect
            float ft = mFocalDistance / viewDir.z;
            Vector3 pFocus = viewDir * ft;
            Vector3 pDxFocus = dxViewDir * ft;
            Vector3 pDyFocus = dyViewDir * ft;
            Vector2 lensSample = mLensRadius *
                uniformSampleDisk(sample.lensU1, sample.lensU2);
            Vector3 viewOrigin(lensSample.x, lensSample.y, 0.0f);
            ray->o = ray->dxOrigin = ray->dyOrigin = 
                mOrientation * viewOrigin + mPosition;
            ray->d = mOrientation * normalize(pFocus - viewOrigin);
            ray->dxDir = mOrientation * normalize(pDxFocus - viewOrigin);
            ray->dyDir = mOrientation * normalize(pDyFocus - viewOrigin);
        }
        ray->mint = 0.0f;
        ray->maxt = INFINITY;
        ray->depth = 0;
        ray->hasDifferential = true;
        return 1.0f;
    }

    OrthographicCamera::OrthographicCamera(const Vector3& position, 
        const Quaternion& orientation, float zn, float zf, 
        float filmWidth, Film* film):
        Camera(position, orientation, zn, zf, film),
        mFilmWidth(filmWidth) {
        mFilmHeight = mFilmWidth / mAspectRatio;
        mProj = matrixOrthoLHD3D(mFilmWidth, mFilmHeight, mZNear, mZFar);
        update();
    }

    float OrthographicCamera::generateRay(const Sample& sample, 
        RayDifferential* ray) const {
        float invXRes = mFilm->getInvXResolution();
        float invYRes = mFilm->getInvYResolution();
        float xNDC = +2.0f * sample.imageX * invXRes - 1.0f;
        float yNDC = -2.0f * sample.imageY * invYRes + 1.0f;
        float dxNDC = +2.0f * (sample.imageX + 1.0f) * invXRes - 1.0f;
        float dyNDC = -2.0f * (sample.imageY + 1.0f) * invYRes + 1.0f;

        // from NDC space to view space
        float xView = 0.5f * mFilmWidth * xNDC;
        float yView = 0.5f * mFilmHeight * yNDC;
        float dxView = 0.5f * mFilmWidth * dxNDC;
        float dyView = 0.5f * mFilmHeight * dyNDC;
        Vector3 pView(xView, yView, 0.0f);
        Vector3 pDxView(dxView, yView, 0.0f);
        Vector3 pDyView(xView, dyView, 0.0f);
        Vector3 viewDir(0.0f, 0.0f, 1.0f);
        // view space to world space
        ray->o = mPosition + mOrientation * pView;
        ray->dxOrigin = mPosition + mOrientation * pDxView;
        ray->dyOrigin = mPosition + mOrientation * pDyView;
        ray->d = ray->dxDir = ray->dyDir = mOrientation * viewDir;
        ray->mint = 0.0f;
        ray->maxt = INFINITY;
        ray->depth = 0;
        ray->hasDifferential = true;
        return 1.0f;
    }


    Camera* PerspectiveCameraCreator::create(const ParamSet& params, 
        Film* film) const {
        Vector3 position = params.getVector3("position");
        Quaternion orientation = getQuaternion(params);
        float fov = params.getFloat("fov", 60.0f);
        float zn = params.getFloat("near_plane", 0.1f);
        float zf = params.getFloat("far_plane", 1000.0f);
        float lensRadius = params.getFloat("lens_radius", 0.0f);
        float focalDistance = params.getFloat("focal_distance", 1e+10f);
        return new PerspectiveCamera(position, orientation, radians(fov), 
            zn, zf, lensRadius, focalDistance, film);
    }


    Camera* OrthographicCameraCreator::create(const ParamSet& params, 
        Film* film) const {
        Vector3 position = params.getVector3("position");
        Quaternion orientation = getQuaternion(params);
        float zn = params.getFloat("near_plane", 0.1f);
        float zf = params.getFloat("far_plane", 1000.0f);
        float filmWidth = params.getFloat("film_width", 35.0f);
        return new OrthographicCamera(position, orientation, zn, zf, 
            filmWidth, film);
    }
}
