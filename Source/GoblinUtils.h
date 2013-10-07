#ifndef GOBLIN_UTILS_H
#define GOBLIN_UTILS_H

#include <cfloat>
#include <cmath>

#ifndef INFINITY
#define INFINITY FLT_MAX
#endif // INFINITY

namespace Goblin {
    const float PI = 3.14159265358979323f;
    const float INV_PI = 0.31830988618379067154f;
    const float INV_TWOPI = 0.15915494309189533577f;

    inline bool isNaN(float f) {
	    // isnan support only in C99/C++11
	    return f != f;
    }

    inline float lerp(float t, float v1, float v2) {
        return (1.0f - t) * v1 + t * v2;
    }

    inline int ceilInt(float f) {
        return (int)ceil(f);
    }

    inline int floorInt(float f) {
        return (int)floor(f);
    }

    inline int roundInt(float f) {
        return floorInt(f + 0.5f);
    }

    inline float radians(float degrees) {
	    return PI * (degrees / 180.0f);
    }

    inline float degrees(float radians) {
	    return (radians / PI) * 180.0f;
    }

    // utils for quadratic equation, solution stored in t1 and t2
    // if it return true
    bool quadratic(float A, float B, float C, float* t1, float* t2);
}

#endif // GOBLIN_UTILS_H