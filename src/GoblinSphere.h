#ifndef GOBLIN_SPHERE_H
#define GOBLIN_SPHERE_H
#include "GoblinGeometry.h"

namespace Goblin {
class Sphere : public Geometry {
public:
    Sphere(float r, size_t numSlices = 30, size_t numStacks = 30);
    ~Sphere() {};
    void init();
    bool intersect(const Ray& ray) const;
    bool intersect(const Ray& ray, float* epsilon, 
        Fragment* fragment) const;
    Vector3 sample(float u1, float u2, Vector3* normal) const;
    Vector3 sample(const Vector3& p, float u1, float u2, 
        Vector3* normal) const;
    float pdf(const Vector3& p, const Vector3& wi) const; 
    float area() const;
    BBox getObjectBound() const;

    size_t getVertexNum() const;
    size_t getFaceNum() const;
    const Vertex* getVertexPtr(size_t index) const;
    const TriangleIndex* getFacePtr(size_t index) const;
private:
    void buildStacks();

private:
    float mRadius;
    size_t mNumSlices;
    size_t mNumStacks;

    VertexList mVertices;
    TriangleList mTriangles;
};

inline float Sphere::area() const {
    return 4.0f * PI * mRadius * mRadius;
}

inline size_t Sphere::getVertexNum() const { 
    return mVertices.size();
}

inline size_t Sphere::getFaceNum() const {
    return mTriangles.size();
}

inline const Vertex* Sphere::getVertexPtr(size_t index) const {
    return &mVertices[index];
}

inline const TriangleIndex* Sphere::getFacePtr(size_t index) const {
    return &mTriangles[index];
}


class ParamSet;
class SceneCache;

Geometry* createSphere(const ParamSet& params, const SceneCache& sceneCache);

}

#endif //GOBLIN_SHPERE_H
