#ifndef GOBLIN_FILTER_H
#define GOBLIN_FILTER_H

#include "GoblinUtils.h"
#include "GoblinParamSet.h"

namespace Goblin {

class Filter {
public:
    Filter(float xWidth, float yWidth);

	virtual ~Filter() = default;

    float getXWidth() const;

    float getYWidth() const;

    virtual float evaluate(float x, float y) const = 0;

    // 2D integration result from -filterWidht->filterWidth
    virtual float getNormalizeTerm() const = 0;

protected:
    float mXWidth, mYWidth;
};


class BoxFilter : public Filter {

public:
    BoxFilter(float xWidth, float yWidth);
    float evaluate(float x, float y) const;
    float getNormalizeTerm() const;
};


class TriangleFilter : public Filter {

public:
    TriangleFilter(float xWidth, float yWidth);
    float evaluate(float x, float y) const;
    float getNormalizeTerm() const;
};


class GaussianFilter : public Filter {
public:
    GaussianFilter(float xWidth, float yWidth, float falloff);

    float evaluate(float x, float y) const;

    float getNormalizeTerm() const;

private:
    float gaussian(float v, float expbase) const;

private:
    const float mAlpha;
    const float mExpX, mExpY;
};


class MitchellFilter: public Filter {
public:
    MitchellFilter(float xWidth, float yWidth, float b, float c) :
		Filter(xWidth, yWidth),
		mInvWidthX(1.0f / xWidth), mInvWidthY(1.0f / yWidth),
		mB(b), mC(c) {}

    float evaluate(float x, float y) const;

    float getNormalizeTerm() const;

private:
    float Mitchell(float x) const;

private:
    const float mInvWidthX, mInvWidthY;
    const float mB, mC;
};

inline Filter* createBoxFilter(const ParamSet& params) {
	Vector2 width = params.getVector2("width", Vector2(1.0f, 1.0f));
	return new BoxFilter(width.x, width.y);
}

inline Filter* createTriangleFilter(const ParamSet& params) {
	Vector2 width = params.getVector2("width", Vector2(1.0f, 1.0f));
	return new TriangleFilter(width.x, width.y);
}

inline Filter* createGaussianFilter(const ParamSet& params) {
	Vector2 width = params.getVector2("width", Vector2(1.0f, 1.0f));
	float falloff = params.getFloat("falloff", 2.0f);
	return new GaussianFilter(width.x, width.y, falloff);
}

inline Filter* createMitchellFilter(const ParamSet& params) {
	Vector2 width = params.getVector2("width", Vector2(1.0f, 1.0f));
	float b = params.getFloat("b", 2.0f);
	float c = params.getFloat("c", 2.0f);
	return new MitchellFilter(width.x, width.y, b, c);
}

}

#endif //GOBLIN_FILTER_H
