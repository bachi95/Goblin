#ifndef GOBLIN_QUATERNION_H
#define GOBLIN_QUATERNION_H

#include "GoblinVector.h"

namespace Goblin {

    class Matrix3;
    class Quaternion {
    public:
        Quaternion();
        Quaternion(float w, float x, float y, float z);
        Quaternion(float w, const Vector3& v);
        Quaternion(const Vector3& axis, float angle);
        Quaternion(const Matrix3& R);
        Quaternion conjugate() const;
        float norm() const;
        Matrix3 toMatrix() const;
        Quaternion operator*(const Quaternion& rhs) const;
        Vector3 operator*(const Vector3& p) const;
        bool operator==(const Quaternion& rhs) const; 
        bool operator!=(const Quaternion& rhs) const;
    public:
        float w;
        Vector3 v;
    }; 

    inline Quaternion::Quaternion() {}

    inline Quaternion::Quaternion(float w, float x, float y, float z):
        w(w), v(x, y, z) {}

    inline Quaternion::Quaternion(float w, const Vector3& v):
        w(w), v(v) {}

    inline Quaternion Quaternion::conjugate() const {
        return Quaternion(w, -v);
    }

    inline float Quaternion::norm() const {
        return w * w + v.squaredLength();
    }

    inline Quaternion Quaternion::operator*(const Quaternion& rhs) const {
        return Quaternion(w * rhs.w - dot(v, rhs.v),
            w * rhs.v + rhs.w * v + cross(v, rhs.v));
    }

    inline bool Quaternion::operator==(const Quaternion& rhs) const {
        return w == rhs.w && v == rhs.v;
    }

    inline bool Quaternion::operator!=(const Quaternion& rhs) const {
        return !operator==(rhs);
    }

    Quaternion normalize(const Quaternion& q);

    std::ostream& operator<<(std::ostream& os, const Quaternion& q);
}

#endif //GOBLIN_QUATERNION_H
