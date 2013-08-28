#include "GoblinCamera.h"
#include "GoblinUtils.h"
#include "GoblinQuaternion.h"
namespace Goblin{

    Camera::Camera() {
        mPosition = Vector3(0.0f, 0.0f, 0.0f);
        mRight = Vector3(1.0f, 0.0f, 0.0f);
        mUp = Vector3(0.0f, 1.0f, 0.0f);
        mLook = Vector3(0.0f, 0.0f, 1.0f);

        mView = Matrix4::Identity;
        mProj = Matrix4::Identity;
    }

    Camera::~Camera() {}

    Matrix4 Camera::view() const { return mView; }

    Matrix4 Camera::proj() const { return mProj; }

    void Camera::setLens(float fovY, float aspect, float zNear, float zFar) {
        mProj = matrixPerspectiveLHD3D(fovY, aspect, zNear, zFar);
    }

    void Camera::strafe(float d) { mPosition += d * mRight; }

    void Camera::walk(float d) { mPosition += d * mLook; }

    void Camera::pitch(float angle) {
        Quaternion R = Quaternion(mRight, angle);
        mUp = R * mUp;
        mLook = R * mLook;
    }

    void Camera::rotateY(float angle) {
        Quaternion R = Quaternion(Vector3::UnitY, angle);
        mRight = R * mRight;
        mUp = R * mUp;
        mLook = R * mLook;
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

}