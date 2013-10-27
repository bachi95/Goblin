#ifndef GOBLIN_RAY_H
#define GOBLIN_RAY_H

#include "GoblinVector.h"
#include "GoblinUtils.h"

namespace Goblin {
    class Ray {
    public:
        Ray();
        Ray(const Vector3& origin, const Vector3& dir, 
            float start, float end = INFINITY, int depth = 0);
        Vector3 operator()(float t) const;
    public:
        Vector3 o;
        Vector3 d;
        mutable float mint, maxt;
        int depth;
    };
    
    inline Ray::Ray(): o(Vector3::Zero), d(Vector3::Zero),
        mint(0), maxt(INFINITY), depth(0) {}

    inline Ray::Ray(const Vector3& origin, const Vector3& dir, 
        float start, float end, int depth): o(origin), d(dir), 
        mint(start), maxt(end), depth(depth) {}

    inline Vector3 Ray::operator()(float t) const {
        return o + t * d;
    }

}

#endif //GOBLIN_RAY_H
