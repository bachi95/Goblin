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

    float BoxFilter::getNormalizeTerm() const {
        return 4.0f * mXWidth * mYWidth;
    }

    TriangleFilter::TriangleFilter(float xWidth, float yWidth):
        Filter(xWidth, yWidth) {}

    float TriangleFilter::evaluate(float x, float y) const {
        return  max(0.0f, mXWidth - fabsf(x)) * max(0.0f, mYWidth - fabsf(y));
    }

    float TriangleFilter::getNormalizeTerm() const {
        // 4 * integrate (integrate (w - x) * (h - y) over 0-w) over 0-h =
        // w * w * h * h
        return mXWidth * mXWidth * mYWidth * mYWidth;
    }

    GaussianFilter::GaussianFilter(float xWidth, float yWidth,
        float falloff): Filter(xWidth, yWidth), mAlpha(falloff),
        mExpX(expf(-falloff * xWidth * xWidth)),
        mExpY(expf(-falloff * yWidth * yWidth)) {}

    float GaussianFilter::evaluate(float x, float y) const {
        return gaussian(x, mExpX) * gaussian(y, mExpY);
    }

    float GaussianFilter::getNormalizeTerm() const {
        // for now just simply use numercal approximation
        // this integration actually involves an erf function
        // that need to be approximated anyway...
        size_t step = 20;
        float deltaX = mXWidth / static_cast<float>(step);
        float deltaY = mYWidth / static_cast<float>(step);
        float result = 0.0f;
        for (size_t i = 0; i < step; ++i) {
            for (size_t j = 0; j < step; ++j) {
                result += 4.0f * deltaX *deltaY *
                    gaussian(i * deltaX, mExpX) *
                    gaussian(j * deltaY, mExpY);
            }
        }
        return result;
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

    float MitchellFilter::getNormalizeTerm() const {
        return 4.0f * ((12 - 9 * mB - 6 * mC) / 4 +
            (-18 + 12 * mB +6 * mC) / 3 + (6 - 2 * mB) +
            15 * (-mB - 6 * mB) / 4 + 7 *(6 * mB + 30 * mC) / 3 +
            3 * (-12 * mB - 48 * mC) / 2 + (8 * mB + 24 * mC)) / 6.0f;
    }

    float MitchellFilter::Mitchell(float x) const {
        x = fabs(2.0f * x);
        float m = 0.0f;
        if (x > 1.0f) {
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
