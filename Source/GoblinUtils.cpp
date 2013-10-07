#include "GoblinUtils.h"
#include <algorithm>

namespace Goblin {

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