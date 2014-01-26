#include "GoblinUtils.h"
#include <ctime>
#include <limits>
#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>


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