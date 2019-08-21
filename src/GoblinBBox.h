#ifndef GOBLIN_BBOX_H
#define GOBLIN_BBOX_H

#include "GoblinVector.h"
#include "GoblinUtils.h"

namespace Goblin {
class Ray;

class BBox {
public:
    BBox() :
		pMin(INFINITY, INFINITY, INFINITY),
		pMax(-INFINITY, -INFINITY, -INFINITY)
	{}

    BBox(const Vector3& p) : pMin(p), pMax(p)
	{}

    BBox(const Vector3& p1, const Vector3& p2) :
		pMin(std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z)),
		pMax(std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z))
	{}

    BBox expand(const Vector3& p);

    BBox expand(const BBox& b);

    BBox expand(float f);

    bool contain(const Vector3& p) const;

    bool intersect(const Ray& ray) const;

    bool intersect(const Ray& ray, float* tMin, float* tMax) const;

    int longestAxis() const;

	Vector3 center() const {
		return 0.5f * (pMin + pMax);
	}

	const Vector3& operator[](int i) const {
		return (&pMin)[i];
	}

	Vector3& operator[](int i) {
		return (&pMin)[i];
	}

	void getBoundingSphere(Vector3* center, float* radius) const {
		*center = 0.5f * (pMin + pMax);
		*radius = length(pMax - pMin);
	}

public:
    Vector3 pMin;
    Vector3 pMax;
};

}

#endif //GOBLIN_BBOX_H