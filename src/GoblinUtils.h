#ifndef GOBLIN_UTILS_H
#define GOBLIN_UTILS_H

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <list>
#include <string>
#include <utility>
#include <vector>
#include <map>
#include <memory>
#include <thread>

#include "GoblinVector.h"

#include <stdint.h> 

#if defined(_WIN32) || defined(_WIN64)
#include <float.h>
#define isinf(f) (!_finite((f)))
#endif

#ifndef INFINITY
#define INFINITY FLT_MAX
#endif // INFINITY

namespace Goblin {

class Color;
class Camera;
class ParamSet;
class Quaternion;
class Renderer;
class Scene;
class Transform;
typedef std::shared_ptr<Camera> CameraPtr;
typedef std::shared_ptr<Renderer> RendererPtr;
typedef std::shared_ptr<Scene> ScenePtr;

const float PI = 3.14159265358979323f;
const float TWO_PI = 6.28318530718f;
const float INV_PI = 0.31830988618379067154f;
const float INV_TWOPI = 0.15915494309189533577f;

inline bool isEqual(float a, float b, const float epsilon = 1e-7f) {
    return fabs(a - b) <= epsilon;
}

inline bool isNaN(float f) {
	// isnan support only in C99/C++11
	return f != f;
}

inline bool isInfinite(const Vector2& v) {
    return v.x >= INFINITY || v.y >= INFINITY;
}

inline bool isInfinite(const Vector3& v) {
    return v.x >= INFINITY || v.y >= INFINITY || v.z >= INFINITY;
}

inline bool isInfinite(const Vector4& v) {
    return v.x >= INFINITY || v.y >= INFINITY || v.z >= INFINITY ||
        v.w >= INFINITY;
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

template <typename T>
T lerp(float t, const T& v1, const T& v2) {
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
    if (root) {
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

inline float sphericalTheta(const Vector3& v) {
    return acos(clamp(v.z, -1.0f, 1.0f));
}

inline float sphericalPhi(const Vector3& v) {
    float phi = atan2(v.y, v.x);
    return phi < 0.0f ? phi + TWO_PI : phi;
}

// solve the 2x2 linear equation B = A * [x, y]
inline bool solve2x2LinearSystem(const float A[2][2], const float B[2],
    float* x, float* y) {
    float det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    if (fabs(det) < 1e-10f) {
        return false;
    }
    *x = (+A[1][1] * B[0] - A[0][1] * B[1]) / det;
    *y = (-A[1][0] * B[0] + A[0][0] * B[1]) / det;
    if (isNaN(*x) || isNaN(*y)) {
        return false;
    }
    return true;
}

// utils that let you form a local coordinate with feed in random axis,
// ex: form a local coordinate based on surface normal vector
void coordinateAxises(const Vector3& a1, Vector3* a2, Vector3* a3);

Quaternion getQuaternion(const ParamSet& params);

Transform getTransform(const ParamSet& params);

// utils for quadratic equation, solution stored in t1 and t2
// if it return true
bool quadratic(float A, float B, float C, float* t1, float* t2);

void drawLine(const Vector2& p0, const Vector2& p1, Color* buffer, 
    int xRes, int yRes, const Color& color); 

void drawPoint(const Vector2& p, Color* buffer, int xRes, int yRes,
    const Color& color, int radius = 1);

// random number generator utils
class RNGImp;
// TODO consider replacing this with PCG
class RNG {
public:
    RNG();
    ~RNG();
    float randomFloat() const;
    uint32_t randomUInt() const;
private:
    RNGImp* mRNGImp;
};

// get first N prime numbers sequence
void getPrimes(size_t N, std::vector<uint32_t>& primes);

inline unsigned int getMaxThreadNum() {
	return std::thread::hardware_concurrency();
}
}

#endif // GOBLIN_UTILS_H
