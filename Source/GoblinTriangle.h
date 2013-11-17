#ifndef GOBLIN_TRIANGLE_H
#define GOBLIN_TRIANGLE_H

#include "GoblinGeometry.h"

namespace Goblin {
    class ObjMesh;
    class Triangle : public Geometry {
    public:
        Triangle(ObjMesh* parentMesh, size_t index);
        bool intersect(const Ray& ray);
        bool intersect(const Ray& ray, float* epsilon, 
            Intersection* intersection);
        BBox getObjectBound();
    private:
        ObjMesh* mParentMesh;
        size_t mIndex;
    };
}

#endif //GOBLIN_TRIANGLE_H