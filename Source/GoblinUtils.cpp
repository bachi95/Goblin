#include "GoblinUtils.h"
#include "GoblinVector.h"

#include <ctime>
#include <limits>

#include <boost/generator_iterator.hpp>
#include <boost/random.hpp>



namespace Goblin {

    // random number generator utils from boost
    typedef boost::mt19937 RNGType;
    const static RNGType rng(static_cast<uint32_t>(time(NULL)));
    const static boost::uniform_real<float> uniRealDist(0.0f, 1.0f);
    const static boost::uniform_int<uint32_t> uniU32Dist(0, 
        numeric_limits<uint32_t>::max() );

    static boost::variate_generator< RNGType, boost::uniform_real<float> >
        randomFloatGen(rng, uniRealDist);

    static boost::variate_generator< RNGType, boost::uniform_int<uint32_t> >
        randomUIntGen(rng, uniU32Dist);

    float randomFloat() {
        return randomFloatGen();
    }

    uint32_t randomUInt() {
        return randomUIntGen();
    }

    void coordinateAxises(const Vector3& a1, Vector3* a2, Vector3* a3) {
        // in case you throw in case like a1 = Vector3(0, 1, 0)
        if(fabsf(a1.x) > fabsf(a1.y)) {
            float invLen = 1.0f / sqrtf(a1.x * a1.x + a1.z * a1.z);
            // dot(a1, a2) = 0 <-> penpenticular to each other
            *a2 = Vector3(-a1.z * invLen, 0.0f, a1.x * invLen);
        } else {
            float invLen = 1.0f / sqrtf(a1.y * a1.y + a1.z * a1.z);
            *a2 = Vector3(0.0f, -a1.z * invLen, a1.y * invLen);
        }
        *a3 = cross(a1, *a2);
    }

    bool quadratic(float A, float B, float C, float* t1, float* t2) {
        float discriminant = B * B - 4.0f * A * C;
        if(discriminant < 0.0f) {
            return false;
        }
        float rootDiscrim = sqrt(discriminant);
        float q;
        // a small trick to avoid numeric error introduced from naive
        // implementation when B or -B close to rootDiscrim
        if(B < 0) {
            q = -0.5f * (B - rootDiscrim);
        } else {
            q = -0.5f * (B + rootDiscrim);
        }
        *t1 = q / A;
        *t2 = C / q;
        if(*t1 > *t2) {
            std::swap(*t1, *t2);
        }
        return true;
    }
}