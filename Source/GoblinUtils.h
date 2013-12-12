#ifndef GOBLIN_UTILS_H
#define GOBLIN_UTILS_H

#include <cfloat>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

// VS2008 doesn't have stdint.h
#if (defined _MSC_VER && _MSC_VER < 1500)
typedef __int8 int8_t; 
typedef unsigned __int8 uint8_t; 
typedef __int16 int16_t; 
typedef unsigned __int16 uint16_t;  
typedef __int32 int32_t; 
typedef unsigned __int32 uint32_t; 
typedef __int64 int64_t; 
typedef unsigned __int64 uint64_t;  
#else 
#include <stdint.h> 
#endif

#ifndef INFINITY
#define INFINITY FLT_MAX
#endif // INFINITY

namespace Goblin {
    using std::min;
    using std::max;
    using std::swap;

    const float PI = 3.14159265358979323f;
    const float TWO_PI = 6.28318530718f;
    const float INV_PI = 0.31830988618379067154f;
    const float INV_TWOPI = 0.15915494309189533577f;

    inline bool isNaN(float f) {
	    // isnan support only in C99/C++11
	    return f != f;
    }

    inline float clamp(float f, float low, float high) {
        if (f < low) {
            return low;
        }
        else if (f > high) {
            return high;
        }
        else return f;
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