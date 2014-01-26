#ifndef GOBLIN_CAMERA_H
#define GOBLIN_CAMERA_H

#include "GoblinMatrix.h"
#include "GoblinVector.h"
#include "GoblinQuaternion.h"

namespace Goblin {
    class Film;
    class Sample;
    class Ray;
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

        Film* getFilm();
        float generateRay(const Sample& sample, Ray* ray);

        bool isUpdated() const;
        void update();

        const Matrix4& getViewMatrix();
        const Matrix4& getWorldMatrix();
        const Matrix4& getProjectionMatrix();

        const Vector3 getLook() const;
        const Vector3 getUp() const;
        const Vector3 getRight() const;

        void roll(float angle);
        void pitch(float angle);
        void yaw(float angle);

        void rotateX(float angel);
        void rotateY(float angle);
        void rotateZ(float angle);

        void rotate(const Vector3& axis, float angle);
        void translate(const Vector3& d);

    private:
        Vector3 mPosition;
        Quaternion mOrientation;
        float mZNear;
        float mZFar;
        float mFOV;
        float mAspectRatio;
        Film* mFilm;

        Matrix4 mWorld;
        Matrix4 mView;
        Matrix4 mProj;

        bool mIsUpdated;
    };

    inline void Camera::setPosition(const Vector3& position) {
        mPosition = position;
        mIsUpdated = false;
    }

    inline const Vector3& Camera::getPosition() const {
        return mPosition;
    }

    inline void Camera::setOrientation(const Quaternion& orientation) {
        mOrientation = orientation;
        mIsUpdated = false;
    }

    inline const Quaternion& Camera::getOrientation() const {
        return mOrientation;
    }

    inline Film* Camera::getFilm(){
        return mFilm;
    }

    inline bool Camera::isUpdated() const {
        return mIsUpdated;
    }
}

#endif //GOBLIN_CAMERA_H
