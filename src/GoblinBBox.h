#ifndef GOBLIN_BBOX_H
#define GOBLIN_BBOX_H

#include "GoblinVector.h"
#include "GoblinUtils.h"

namespace Goblin {
class Ray;

class BBox {
public:
    BBox();
    BBox(const Vector3& p);
    BBox(const Vector3& p1, const Vector3& p2);
    BBox expand(const Vector3& p);
    BBox expand(const BBox& b);
    BBox expand(float f);
    bool contain(const Vector3& p) const;
    bool intersect(const Ray& ray) const;
    bool intersect(const Ray& ray, float* tMin, float* tMax) const;
    int longestAxis() const;
    Vector3 center() const;
    const Vector3& operator[](int i) const;
    Vector3& operator[](int i);
    void getBoundingSphere(Vector3* center, float* radius) const;

public:
    Vector3 pMin;
    Vector3 pMax;
};

inline BBox::BBox():
    pMin(INFINITY, INFINITY, INFINITY),
    pMax(-INFINITY, -INFINITY, -INFINITY) {}

inline BBox::BBox(const Vector3& p): pMin(p), pMax(p) {}

inline BBox::BBox(const Vector3& p1, const Vector3& p2):
    pMin(min(p1.x, p2.x), min(p1.y, p2.y), min(p1.z, p2.z)),
    pMax(max(p1.x, p2.x), max(p1.y, p2.y), max(p1.z, p2.z)) {}

inline Vector3 BBox::center() const {
    return 0.5f * (pMin + pMax);
}

inline const Vector3& BBox::operator[](int i) const {
    assert(i == 0 || i == 1);
    return (&pMin)[i];
}

inline Vector3& BBox::operator[](int i) {
    assert(i == 0 || i == 1);
    return (&pMin)[i];
}

inline void BBox::getBoundingSphere(Vector3* center, float* radius) const {
    *center = 0.5f * (pMin + pMax);
    *radius = length(pMax - pMin);
}

}

#endif //GOBLIN_BBOX_H
