#ifndef GOBLIN_DISK_H
#define GOBLIN_DISK_H
#include "GoblinGeometry.h"

namespace Goblin {
class Disk : public Geometry {
public:
    Disk(float r, size_t numSlices = 30);
    ~Disk() {};
    void init();
    bool intersect(const Ray& ray) const;
    bool intersect(const Ray& ray, float* epsilon,
        Fragment* fragment) const;
    Vector3 sample(float u1, float u2, Vector3* normal) const;
    float area() const;
    BBox getObjectBound() const;

    size_t getVertexNum() const;
    size_t getFaceNum() const;
    const Vertex* getVertexPtr(size_t index) const;
    const TriangleIndex* getFacePtr(size_t index) const;
private:
    void buildSlices();

private:
    float mRadius;
    size_t mNumSlices;

    VertexList mVertices;
    TriangleList mTriangles;
};

inline float Disk::area() const {
    return PI * mRadius * mRadius;
}

inline size_t Disk::getVertexNum() const {
    return mVertices.size();
}

inline size_t Disk::getFaceNum() const {
    return mTriangles.size();
}

inline const Vertex* Disk::getVertexPtr(size_t index) const {
    return &mVertices[index];
}

inline const TriangleIndex* Disk::getFacePtr(size_t index) const {
    return &mTriangles[index];
}

class ParamSet;
class SceneCache;

Geometry* createDisk(const ParamSet& params, const SceneCache& sceneCache);

}

#endif //GOBLIN_DISK_H