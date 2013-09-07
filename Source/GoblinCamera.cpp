#include "GoblinCamera.h"
#include "GoblinFilm.h"
#include "GoblinRay.h"
#include "GoblinSampler.h"
#include "GoblinUtils.h"

namespace Goblin{

    Camera::Camera(const Vector3& position, const Quaternion& orientation,
        float fov, float zn, float zf, Film* film):
        mPosition(position), mOrientation(orientation),
        mFOV(fov), mZNear(zn), mZFar(zf), mFilm(film) {
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

        mRight = Vector3(1.0f, 0.0f, 0.0f);
        mUp = Vector3(0.0f, 1.0f, 0.0f);
        mLook = Vector3(0.0f, 0.0f, 1.0f);

        mProj = matrixPerspectiveLHD3D(PI * 0.25f, 1.3333f, 0.0001f, 1000.0f);
        mFilm = NULL;
        update();
    }

    Camera::~Camera() {
        if(mFilm != NULL) {
            delete mFilm;
            mFilm = NULL;
        }
    }

    float Camera::generateRay(const CameraSample& sample, Ray* ray) {
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
        ray->o = Vector3::Zero;
        ray->d = normalize(Vector3(xView, yView, zView));
        ray->mint = 0.0f;
        ray->maxt = INFINITY;
        ray->depth = 0;

        //TODO transform the view space ray to world space
        ray->o = mPosition;
        ray->d = mOrientation * ray->d; 
        return 1.0f;
    }

    Matrix4 Camera::view() const { 
        return mView; 
    }

    Matrix4 Camera::proj() const { 
        return mProj; 
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
    

    void Camera::setLens(float fovY, float aspect, float zNear, float zFar) {
        mProj = matrixPerspectiveLHD3D(fovY, aspect, zNear, zFar);
    }

    void Camera::strafe(float d) { 
        mPosition += d * mRight; 
        mIsUpdated = false;
    }

    void Camera::walk(float d) { 
        mPosition += d * mLook; 
        mIsUpdated = false;
    }

    void Camera::pitch(float angle) {
        Quaternion R = Quaternion(mRight, angle);
        mUp = R * mUp;
        mLook = R * mLook;

        Vector3 right = mOrientation * Vector3::UnitX;
        rotate(right, angle);
        mIsUpdated = false;
    }

    void Camera::rotateY(float angle) {
        Quaternion R = Quaternion(Vector3::UnitY, angle);
        mRight = R * mRight;
        mUp = R * mUp;
        mLook = R * mLook;
        rotate(Vector3::UnitY, angle);
        mIsUpdated = false;
    }
    
    void Camera::rotate(const Vector3& axis,float angle) {
        mOrientation = normalize(Quaternion(axis, angle) * mOrientation);
        mIsUpdated = false;
    } 

    void Camera::rebuildView() {
        mLook.normalize();

        mUp = normalize(cross(mLook, mRight));
        mRight = normalize(cross(mUp, mLook));

        float x = -dot(mPosition, mRight);
        float y = -dot(mPosition, mUp);
        float z = -dot(mPosition, mLook);

        mView = Matrix4(
            mRight.x, mRight.y, mRight.z,    x,
               mUp.x,    mUp.y,    mUp.z,    y,
             mLook.x,  mLook.y,  mLook.z,    z,
                0.0f,     0.0f,     0.0f, 1.0f);
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

}
