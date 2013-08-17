#ifndef GOBLIN_RAY_H
#define GOBLIN_RAY_H

#include "GoblinVector.h"
#include "GoblinUtils.h"

namespace Goblin {
    class Ray {
    public:
        Ray();
        Ray(const Vector3& origin, const Vector3& dir, 
            float start, float end = INFINITY, float depth = 0);
    public:
        Vector3 o;
        Vector3 d;
        float mint, maxt;
        float depth;
    };
    
    inline Ray::Ray(): o(Vector3::Zero), d(Vector3::Zero),
        mint(0), maxt(INFINITY), depth(0) {}

    inline Ray::Ray(const Vector3& origin, const Vector3& dir, 
        float start, float end, float depth): o(origin), d(dir), 
        mint(start), maxt(end), depth(depth) {}

}

#endif //GOBLIN_RAY_H
