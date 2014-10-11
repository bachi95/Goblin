#ifndef GOBLIN_SPHERE_H
#define GOBLIN_SPHERE_H
#include "GoblinGeometry.h"
#include "GoblinFactory.h"

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
        float area() const;
        BBox getObjectBound();

    private:
        void buildStacks();

    private:
        float mRadius;
        size_t mNumSlices;
        size_t mNumStacks;
    };

    inline float Sphere::area() const {
        return 4.0f * PI * mRadius * mRadius;
    }


    class ParamSet;
    class SceneCache;

    class SphereGeometryCreator : 
        public Creator<Geometry, const ParamSet&, const SceneCache&> {
    public:
        Geometry* create(const ParamSet& params, 
            const SceneCache& sceneCache) const;
    };
}

#endif //GOBLIN_SHPERE_H
