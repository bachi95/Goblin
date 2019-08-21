#ifndef GOBLIN_TRANSFORM_H
#define GOBLIN_TRANSFORM_H

#include "GoblinVector.h"
#include "GoblinMatrix.h"
#include "GoblinQuaternion.h"

namespace Goblin {
class Ray;
class BBox;

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
    const Matrix4& getMatrix() const;
    const Matrix4& getInverse() const;

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
    BBox onBBox(const BBox& b) const;

    Vector3 invertPoint(const Vector3& p) const;
    Vector3 invertNormal(const Vector3& n) const;
    Vector3 invertVector(const Vector3& v) const;
    Ray invertRay(const Ray& ray) const;
    BBox invertBBox(const BBox& b) const;

    bool isUpdated() const;
    void update() const;

private:
    mutable Matrix4 mCachedMatrix;
    mutable Matrix4 mCachedInverse;
    Vector3 mPosition;
    float mPad;
    Quaternion mOrientation;
    Vector3 mScale;
    mutable bool mIsUpdated;
};
}

#endif //GOBLIN_TRANSFORM_H