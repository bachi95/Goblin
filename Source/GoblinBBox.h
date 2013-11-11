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
        bool intersect(const Ray& ray);

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
}

#endif //GOBLIN_BBOX_H
