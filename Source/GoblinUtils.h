#ifndef GOBLIN_UTILS_H
#define GOBLIN_UTILS_H

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <map>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

using namespace std;
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
    using std::vector;
    using std::string;
    using std::map;
    using std::min;
    using std::max;
    using std::swap;
    using std::cout;
    using std::cerr;
    using std::endl;
    using std::pair;

    class Color;
    class Vector2;
    class Vector3;
    class Camera;
    class Scene;
    typedef boost::shared_ptr<Camera> CameraPtr;
    typedef boost::shared_ptr<Scene> ScenePtr;

    const float PI = 3.14159265358979323f;
    const float TWO_PI = 6.28318530718f;
    const float INV_PI = 0.31830988618379067154f;
    const float INV_TWOPI = 0.15915494309189533577f;

    inline bool isNaN(float f) {
	    // isnan support only in C99/C++11
	    return f != f;
    }

    inline bool isPowerOf2(uint32_t n) {
        return (n & (n - 1)) == 0;
    }

    inline uint32_t roundUpPow2(uint32_t n) {
        n--;
        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        return n + 1;
    }

    inline float log2(float n) {
        static float invLog2 = 1.0f / logf(2.0f);
        return logf(n) * invLog2;
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

    inline int clamp(int f, int low, int high) {
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

    inline int roundToSquare(int n, int* root = NULL) {
        int s = ceilInt(sqrt(static_cast<float>(n)));
        if(root) {
            *root = s;
        }
        return s * s;
    }

    inline float radians(float degrees) {
	    return PI * (degrees / 180.0f);
    }

    inline float degrees(float radians) {
	    return (radians / PI) * 180.0f;
    }

    // utils that let you form a local coordinate with feed in random axis,
    // ex: form a local coordinate based on surface normal vector
    void coordinateAxises(const Vector3& a1, Vector3* a2, Vector3* a3);

    // utils for quadratic equation, solution stored in t1 and t2
    // if it return true
    bool quadratic(float A, float B, float C, float* t1, float* t2);

    void drawLine(const Vector2& p0, const Vector2& p1, Color* buffer, 
        int xRes, int yRes, const Color& color); 

    void drawPoint(const Vector2& p, Color* buffer, int xRes, int yRes,
        const Color& color, int radius = 1);

    // random number generator utils from boost
    class RNGImp;

    class RNG {
    public:
        RNG();
        ~RNG();
        float randomFloat() const;
        uint32_t randomUInt() const;
    private:
        RNGImp* mRNGImp;
    };

    class NullType {};
}

#endif // GOBLIN_UTILS_H
