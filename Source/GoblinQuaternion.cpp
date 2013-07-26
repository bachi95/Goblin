#include "GoblinQuaternion.h"
#include <cmath>

namespace Goblin {
    Quaternion::Quaternion(const Vector3& axis, float radians) {
        float t = radians * 0.5f;
        Vector3 unitAxis = normalize(axis);
        float sinT = sin(t);
        w = cos(t);
        v = unitAxis * sinT;
    } 

    // the NVIDIA optimized version of qpq^-1
    Vector3 Quaternion::operator*(const Vector3& p) const {
        Vector3 uv = cross(v, p);
        Vector3 uuv = cross(v, uv);
        uv *= (2.0f * w);
        uuv *= 2.0f;

        return p + uv + uuv;
    }

    Quaternion normalize(const Quaternion& q) {
        float inv = 1.0f / sqrt(q.norm());
        Quaternion result(q);
        result.v *= inv;
        result.w *= inv;
        return result;
    }

    std::ostream& operator<<(std::ostream& os, const Quaternion& q) {
        os << "Quaternion(" << q.w << ", " << 
            q.v.x << ", " << q.v.y << ", " << q.v.z <<")";
        return os;
    }
}
