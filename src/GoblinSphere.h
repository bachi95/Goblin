#ifndef GOBLIN_SPHERE_H
#define GOBLIN_SPHERE_H
#include "GoblinGeometry.h"

namespace Goblin {

class Sphere : public Geometry {
public:
    Sphere(float r);

	~Sphere() = default;

    bool intersect(const Ray& ray, float* epsilon, 
        Fragment* fragment) const override;

	bool occluded(const Ray& ray) const override;

    Vector3 sample(float u1, float u2, Vector3* normal) const override;

    Vector3 sample(const Vector3& p, float u1, float u2, 
        Vector3* normal) const override;

    float pdf(const Vector3& p, const Vector3& wi) const override;

    float area() const override {
		return 4.0f * PI * mRadius * mRadius;
	}

    BBox getObjectBound() const override;

private:
    float mRadius;
};

class ParamSet;
class SceneCache;

Geometry* createSphere(const ParamSet& params, const SceneCache& sceneCache);

}

#endif //GOBLIN_SHPERE_H