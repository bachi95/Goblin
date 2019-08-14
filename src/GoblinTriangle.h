#ifndef GOBLIN_TRIANGLE_H
#define GOBLIN_TRIANGLE_H

#include "GoblinGeometry.h"

namespace Goblin {
class ObjMesh;
class Triangle : public Geometry {
public:
    Triangle(const ObjMesh* parentMesh): mParentMesh(parentMesh), mIndex(0) {}

    Triangle(const ObjMesh* parentMesh, size_t index);

    void setIndex(size_t index);

    bool intersect(const Ray& ray, float* epsilon,
        Fragment* fragment) const override;

	bool occluded(const Ray& ray) const override;

    Vector3 sample(float u1, float u2, Vector3* normal) const override;

    float area() const override;

    BBox getObjectBound() const override;

private:
    const ObjMesh* mParentMesh;
    size_t mIndex;
};

inline void Triangle::setIndex(size_t index) { mIndex = index; }

}

#endif //GOBLIN_TRIANGLE_H
