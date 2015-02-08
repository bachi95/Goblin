#ifndef GOBLIN_CAMERA_H
#define GOBLIN_CAMERA_H

#include "GoblinMatrix.h"
#include "GoblinVector.h"
#include "GoblinQuaternion.h"
#include "GoblinFactory.h"

namespace Goblin {
    class Film;
    class Sample;
    class Ray;
    class Camera {
    public:
        Camera(const Vector3& position, const Quaternion& orientation,
            float zn, float zf, Film* film);
        virtual ~Camera();

        void setPosition(const Vector3& position);
        const Vector3& getPosition() const;
        void setOrientation(const Quaternion& orientation);
        const Quaternion& getOrientation() const;

        Film* getFilm();
        virtual float generateRay(const Sample& sample, Ray* ray) const = 0;

        bool isUpdated() const;
        void update();

        const Matrix4& getWorldMatrix();
        const Matrix4& getViewMatrix();
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

        Vector3 worldToScreen(const Vector3& pWorld) const;

    protected:
        Vector3 mPosition;
        Quaternion mOrientation;
        float mZNear;
        float mZFar;
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

    inline const Matrix4& Camera::getProjectionMatrix() {
        return mProj;
    }

    inline const Vector3 Camera::getLook() const {
        return mOrientation * Vector3::UnitZ;
    }

    inline const Vector3 Camera::getUp() const {
        return mOrientation * Vector3::UnitY;
    }

    inline const Vector3 Camera::getRight() const {
        return mOrientation * Vector3::UnitX;
    }

    inline void Camera::roll(float angle) {
        rotate(mOrientation * Vector3::UnitZ, angle);
    }

    inline void Camera::pitch(float angle) {
        rotate(mOrientation * Vector3::UnitX, angle);
    }

    inline void Camera::yaw(float angle) {
        rotate(mOrientation * Vector3::UnitY, angle);
    }

    inline void Camera::rotateX(float angle) {
        rotate(Vector3::UnitX, angle);
    }

    inline void Camera::rotateY(float angle) {
        rotate(Vector3::UnitY, angle);
    }

    inline void Camera::rotateZ(float angle) {
        rotate(Vector3::UnitZ, angle);
    }
    
    inline void Camera::rotate(const Vector3& axis,float angle) {
        mOrientation = normalize(Quaternion(axis, angle) * mOrientation);
        mIsUpdated = false;
    } 

    inline void Camera::translate(const Vector3& d) {
        mPosition += d;
        mIsUpdated = false;
    }


    class PerspectiveCamera : public Camera {
    public:
        PerspectiveCamera(const Vector3& position, 
            const Quaternion& orientation,
            float fov, float zn, float zf, 
            float lensRadius, float focalDistance, 
            Film* film);

        float generateRay(const Sample& sample, Ray* ray) const;
    private:
        float mFOV;
        float mLensRadius;
        float mFocalDistance;
    };


    class OrthographicCamera : public Camera {
    public:
        OrthographicCamera(const Vector3& position, 
            const Quaternion& orientation,
            float zn, float zf, float filmWidth, Film* film);

        float generateRay(const Sample& sample, Ray* ray) const;
    private:
        float mFilmWidth;
        float mFilmHeight;
    };


    class PerspectiveCameraCreator : 
        public Creator<Camera, const ParamSet&, Film*> {
    public:
        Camera* create(const ParamSet& params, Film* film) const;
    };


    class OrthographicCameraCreator : 
        public Creator<Camera, const ParamSet&, Film*> {
    public:
        Camera* create(const ParamSet& params, Film* film) const;
    };
}

#endif //GOBLIN_CAMERA_H
