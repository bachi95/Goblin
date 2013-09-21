#ifndef GOBLIN_SPHERE_H
#define GOBLIN_SPHERE_H
#include "GoblinGeometry.h"
#include "GoblinVector.h"

namespace Goblin {
    class Sphere : public Geometry {
    public:
        Sphere(const Vector3& o,float r, 
            size_t numSlices = 30, 
            size_t numStacks = 30);
        ~Sphere() {};
        void init();
    private:
        void buildStacks();

    private:
        Vector3 mOrigin;
        float mRadius;
        size_t mNumSlices;
        size_t mNumStacks;
    };
}

#endif //GOBLIN_SHPERE_H
