#ifndef GOBLIN_TRIANGLE_H
#define GOBLIN_TRIANGLE_H

#include "GoblinGeometry.h"

namespace Goblin {
    class ObjMesh;
    class Triangle : public Geometry {
    public:
        Triangle(const ObjMesh* parentMesh): mParentMesh(parentMesh) {}
        Triangle(const ObjMesh* parentMesh, size_t index);
        void setIndex(size_t index);
        bool intersect(const Ray& ray) const;
        bool intersect(const Ray& ray, float* epsilon, 
            Fragment* fragment) const;
        size_t getVertexNum() const;
        size_t getFaceNum() const;
        const Vertex* getVertexPtr(size_t index) const;
        const TriangleIndex* getFacePtr(size_t index) const;

        Vector3 sample(float u1, float u2, Vector3* normal) const;
        float area() const;
        BBox getObjectBound() const;
    private:
        const ObjMesh* mParentMesh;
        size_t mIndex;
    };

    inline void Triangle::setIndex(size_t index) { mIndex = index; }

    inline size_t Triangle::getVertexNum() const { return 3; }

    inline size_t Triangle::getFaceNum() const { return 1; }
}

#endif //GOBLIN_TRIANGLE_H
