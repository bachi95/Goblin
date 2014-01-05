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

    inline Filter::Filter(float xWidth, float yWidth):
        mXWidth(xWidth), mYWidth(yWidth) {}

    inline float Filter::getXWidth() const {
        return mXWidth;
    }

    inline float Filter::getYWidth() const {
        return mYWidth;
    }


    class BoxFilter : public Filter {
    public:
        BoxFilter(float xWidth, float yWidth);
        float evaluate(float x, float y) const;
    };

    inline BoxFilter::BoxFilter(float xWidth, float yWidth):
        Filter(xWidth, yWidth) {}

    inline float BoxFilter::evaluate(float x, float y) const {
        return 1.0f;
    }


    class TriangleFilter : public Filter {
    public:
        TriangleFilter(float xWidth, float yWidth);
        float evaluate(float x, float y) const;
    };

    inline TriangleFilter::TriangleFilter(float xWidth, float yWidth):
        Filter(xWidth, yWidth) {}

    inline float TriangleFilter::evaluate(float x, float y) const {
        return  max(0.0f, mXWidth - fabs(x)) * max(0.0f, mYWidth - fabs(y));
    }
    
}

#endif //GOBLIN_FILTER_H
