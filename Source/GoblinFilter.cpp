#include "GoblinFilter.h"

namespace Goblin {
    Filter::Filter(float xWidth, float yWidth):
        mXWidth(xWidth), mYWidth(yWidth) {}

    float Filter::getXWidth() const {
        return mXWidth;
    }

    float Filter::getYWidth() const {
        return mYWidth;
    }

    BoxFilter::BoxFilter(float xWidth, float yWidth):
        Filter(xWidth, yWidth) {}

    float BoxFilter::evaluate(float x, float y) const {
        return 1.0f;
    }

    TriangleFilter::TriangleFilter(float xWidth, float yWidth):
        Filter(xWidth, yWidth) {}

    float TriangleFilter::evaluate(float x, float y) const {
        return  max(0.0f, mXWidth - fabsf(x)) * max(0.0f, mYWidth - fabsf(y));
    }

    GaussianFilter::GaussianFilter(float xWidth, float yWidth,
        float falloff): Filter(xWidth, yWidth), mAlpha(falloff),
        mExpX(expf(-falloff * xWidth * xWidth)),
        mExpY(expf(-falloff * yWidth * yWidth)) {}

    float GaussianFilter::evaluate(float x, float y) const {
        return gaussian(x, mExpX) * gaussian(y, mExpY);
    }

    float GaussianFilter::gaussian(float v, float expbase) const {
        return max(0.0f, expf(-mAlpha * v * v) - expbase);
    }

    MitchellFilter::MitchellFilter(float xWidth, float yWidth,
        float b, float c): Filter(xWidth, yWidth), 
        mInvWidthX(1.0f / xWidth), mInvWidthY(1.0f / yWidth),
        mB(b), mC(c) {}

    float MitchellFilter::evaluate(float x, float y) const {
        return Mitchell(x * mInvWidthX) * Mitchell(y * mInvWidthY);
    }

    float MitchellFilter::Mitchell(float x) const {
        x = fabs(2.0f * x);
        float m = 0.0f;
        if(x > 1.0f) {
            m = ((-mB - 6 * mC) * x * x * x + (6 * mB + 30 * mC) * x * x +
                (-12 * mB - 48 * mC) * x + (8 * mB + 24 * mC)) / 6.0f;
        } else {
            m = ((12 - 9 * mB - 6 * mC) * x * x * x + 
                (-18 + 12 * mB + 6 * mC) * x * x +
                (6 - 2 * mB)) / 6.0f;
        }
        return m;
    }
}
