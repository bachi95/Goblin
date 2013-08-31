#ifndef GOBLIN_CAMERA_H
#define GOBLIN_CAMERA_H

#include "GoblinMatrix.h"
#include "GoblinVector.h"
#include "GoblinQuaternion.h"

namespace Goblin {
    class Film;
    class Camera {
    public:
        Camera();
        Camera(const Vector3& position, const Quaternion& orientation,
            float fov, float zn, float zf, Film* film);
        ~Camera();

        void setPosition(const Vector3& position);
        const Vector3& getPosition() const;
        void setOrientation(const Quaternion& orientation);
        const Quaternion& getOrientation() const;

        Matrix4 view() const;
        Matrix4 proj() const;

        void setLens(float fovY, float aspect, float nearm, float far);
        void strafe(float d);
        void walk(float d);
        void pitch(float angle);
        void rotateY(float angle);

        void rebuildView();
    private:
        Vector3 mPosition;
        Vector3 mRight;
        Vector3 mUp;
        Vector3 mLook;

        float mZNear;
        float mZFar;
        float mFOV;
        float mAspectRatio;
        Quaternion mOrientation;

        Matrix4 mView;
        Matrix4 mProj;
        Film* mFilm;
    };

    inline void Camera::setPosition(const Vector3& position) {
        mPosition = position;
    }

    inline const Vector3& Camera::getPosition() const {
        return mPosition;
    }

    inline void Camera::setOrientation(const Quaternion& orientation) {
        mOrientation = orientation;
    }

    inline const Quaternion& Camera::getOrientation() const {
        return mOrientation;
    }
}

#endif //GOBLIN_CAMERA_H
