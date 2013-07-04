#ifndef GOBLIN_QUATERNION_H
#define GOBLIN_QUATERNION_H

#include "GoblinVector.h"

namespace Goblin {

    class Quaternion {
    public:
        Quaternion() {}
        Quaternion(float w, float x, float y, float z) :
            w(w), v(x, y, z) {}
        Quaternion(float w, const Vector3& v) :
            w(w), v(v) {}
        Quaternion(const Vector3& axis, float radians);

        Quaternion conjugate() {
            return Quaternion(w, -v.x, -v.y, -v.z);
        }

        Quaternion operator*(const Quaternion& rhs) const {
            return Quaternion(w * rhs.w - dot(v, rhs.v),
                w * rhs.v + rhs.w * v + cross(v, rhs.v));
        }

        bool operator==(const Quaternion& rhs) const {
            return w == rhs.w && v == rhs.v;
        } 

        bool operator!=(const Quaternion& rhs) const {
            return !operator==(rhs);
        }

    public:
        float w;
        Vector3 v;
    }; 
}

#endif //GOBLIN_QUATERNION_H
