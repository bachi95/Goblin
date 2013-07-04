#include "GoblinQuaternion.h"
#include "GoblinVector.h"
#include <cmath>

namespace Goblin {
    Quaternion::Quaternion(const Vector3& axis, float radians) {
        float t = radians * 0.5f;
        Vector3 unitAxis = normalize(axis);
        float sinT = sin(t);
        w = cos(t);
        v = unitAxis * sinT;
    } 
}
