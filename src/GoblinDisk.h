#ifndef GOBLIN_DISK_H
#define GOBLIN_DISK_H
#include "GoblinGeometry.h"

namespace Goblin {

class Disk : public Geometry {
public:
    Disk(float r);

	~Disk() = default;

    bool intersect(const Ray& ray, float* epsilon,
		Fragment* fragment) const override;

	bool occluded(const Ray& ray) const override;

    Vector3 sample(float u1, float u2, Vector3* normal) const override;

    float area() const override {
		return PI * mRadius * mRadius;
	}

    BBox getObjectBound() const override;

private:
    float mRadius;
};

class ParamSet;
class SceneCache;

Geometry* createDisk(const ParamSet& params, const SceneCache& sceneCache);

}

#endif //GOBLIN_DISK_H