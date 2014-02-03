#ifndef GOBLIN_SPHERE_H
#define GOBLIN_SPHERE_H
#include "GoblinGeometry.h"

namespace Goblin {
    class Sphere : public Geometry {
    public:
        Sphere(float r, size_t numSlices = 30, size_t numStacks = 30);
        ~Sphere() {};
        void init();
        bool intersect(const Ray& ray);
        bool intersect(const Ray& ray, float* epsilon, 
            Fragment* fragment);
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
}

#endif //GOBLIN_SHPERE_H
