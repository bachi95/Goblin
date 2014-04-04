#ifndef GOBLIN_FILTER_H
#define GOBLIN_FILTER_H

#include "GoblinUtils.h"

namespace Goblin {
    class Filter {
    public:
        Filter(float xWidth, float yWidth);
        float getXWidth() const;
        float getYWidth() const;
        virtual float evaluate(float x, float y) const = 0;
    protected:
        float mXWidth, mYWidth;
    };


    class BoxFilter : public Filter {
    public:
        BoxFilter(float xWidth, float yWidth);
        float evaluate(float x, float y) const;
    };


    class TriangleFilter : public Filter {
    public:
        TriangleFilter(float xWidth, float yWidth);
        float evaluate(float x, float y) const;
    };


    class GaussianFilter : public Filter {
    public:
        GaussianFilter(float xWidth, float yWidth, float falloff);
        float evaluate(float x, float y) const;
    private:
        float gaussian(float v, float expbase) const;
    private:
        const float mAlpha;
        const float mExpX, mExpY;
    };


    class MitchellFilter: public Filter {
    public:
        MitchellFilter(float xWidth, float yWidth, float b, float c);
    private:
        float Mitchell(float x) const;
        float evaluate(float x, float y) const;
    private:
        const float mInvWidthX, mInvWidthY;
        const float mB, mC;
    };

}

#endif //GOBLIN_FILTER_H
