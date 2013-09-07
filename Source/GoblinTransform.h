#ifndef GOBLIN_TRANSFORM_H
#define GOBLIN_TRANSFORM_H

#include "GoblinVector.h"
#include "GoblinMatrix.h"
#include "GoblinQuaternion.h"

namespace Goblin {
    class Ray;
    class Transform {
    public:
        Transform();
        Transform(const Vector3& position, const Quaternion& orientation, 
            const Vector3& scale);

        void setPosition(const Vector3& position);
        void setOrientation(const Quaternion& orientation);
        void setScale(const Vector3& scale);
        const Vector3& getPosition() const;
        const Quaternion& getOrientation() const;
        const Vector3& getScale() const;
        const Matrix4& getMatrix();
        const Matrix4& getInverse();

        void roll(float angle);
        void pitch(float angle);
        void yaw(float angle);

        void rotateX(float angle);
        void rotateY(float angle);
        void rotateZ(float angle);

        void rotate(const Vector3& axis, float angle);
        void translate(const Vector3& d);

        Vector3 onPoint(const Vector3& p) const;
        Vector3 onNormal(const Vector3& n) const;
        Vector3 onVector(const Vector3& v) const;
        Ray onRay(const Ray& ray) const;

        bool isUpdated() const;
        void update();

    private:
        Matrix4 mCachedMatrix;
        Matrix4 mCachedInverse;

        Quaternion mOrientation;
        Vector3 mPosition;
        Vector3 mScale;
        bool mIsUpdated;
    };
}

#endif //GOBLIN_TRANSFORM_H