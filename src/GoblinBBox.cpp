#include "GoblinBBox.h"
#include "GoblinRay.h"

namespace Goblin {
BBox BBox::expand(const Vector3& p) {
    pMin = Vector3(min(pMin.x, p.x), min(pMin.y, p.y), min(pMin.z, p.z));
    pMax = Vector3(max(pMax.x, p.x), max(pMax.y, p.y), max(pMax.z, p.z));
    return *this;
}

BBox BBox::expand(const BBox& b) {
    pMin = Vector3(
        min(pMin.x, b.pMin.x),
        min(pMin.y, b.pMin.y),
        min(pMin.z, b.pMin.z));
    pMax = Vector3(
        max(pMax.x, b.pMax.x),
        max(pMax.y, b.pMax.y),
        max(pMax.z, b.pMax.z));
    return *this;
}

BBox BBox::expand(float f) {
    Vector3 delta(f, f, f);
    pMin -= delta;
    pMax += delta;
    return *this;
}

bool BBox::contain(const Vector3& p) const {
    return pMin.x <= p.x && p.x <= pMax.x &&
        pMin.y <= p.y && p.y <= pMax.y &&
        pMin.z <= p.z && p.z <= pMax.z;
}

bool BBox::intersect(const Ray& ray) const {
    float t0 = ray.mint;
    float t1 = ray.maxt;
    for (int i = 0; i < 3; ++i) {
        float invDir = 1.0f / ray.d[i];
        float tNear = (pMin[i] - ray.o[i]) * invDir;
        float tFar = (pMax[i] - ray.o[i]) * invDir;
        // ray might travel in/away
        if (tNear > tFar) {
            swap(tNear, tFar);
        }
        t0 = (tNear > t0) ? tNear : t0;
        t1 = (tFar < t1) ? tFar : t1;
        if (t0 > t1) {
            return false;
        }
    }
    return true;
}

bool BBox::intersect(const Ray& ray, float* tMin, float* tMax) const {
    float t0 = ray.mint;
    float t1 = ray.maxt;
    for (int i = 0; i < 3; ++i) {
        float invDir = 1.0f / ray.d[i];
        float tNear = (pMin[i] - ray.o[i]) * invDir;
        float tFar = (pMax[i] - ray.o[i]) * invDir;
        // ray might travel in/away
        if (tNear > tFar) {
            swap(tNear, tFar);
        }
        t0 = (tNear > t0) ? tNear : t0;
        t1 = (tFar < t1) ? tFar : t1;
        if (t0 > t1) {
            return false;
        }
    }
    *tMin = t0;
    *tMax = t1;
    return true;
}

int BBox::longestAxis() const {
    Vector3 delta = pMax - pMin;
    if (delta.x > delta.y && delta.x > delta.z) {
        return 0;
    } else if (delta.y > delta.z) {
        return 1;
    } else {
        return 2;
    }
}

} // namespace Goblin