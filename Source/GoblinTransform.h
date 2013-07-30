#ifndef GOBLIN_TRANSFORM_H
#define GOBLIN_TRANSFORM_H

#include "GoblinVector.h"
#include "GoblinMatrix.h"
#include "GoblinQuaternion.h"

namespace Goblin {

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

        Vector3 onPoint(const Vector3& p) const;
        Vector3 onNormal(const Vector3& n) const;
        Vector3 onVector(const Vector3& v) const;
        //Ray onRay(const Ray& ray) const;

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