#ifndef GOBLIN_PRIMITIVE_H
#define GOBLIN_PRIMITIVE_H

#include "GoblinTransform.h"
#include "GoblinGeometry.h"
#include "GoblinMaterial.h"
#include "GoblinBBox.h"
#include "GoblinUtils.h"
#include "GoblinLight.h"

#include <vector>
#include <exception>

namespace Goblin {

class Ray;
class RayDifferential;
/* temp notes:
model = geometry + material
instance = transform + any kind of primitive

they all implement:
intersect
model::intersect -> return geometry.intersect + material
instance::intersect -> transform the ray to object space
	then mPrimitive.intersect
*/

class Primitive;
typedef std::vector<const Primitive*> PrimitiveList;

struct Intersection {
	Intersection() : primitive(nullptr) {}

	Color Le(const Vector3& outDirection);

	const Light* getLight() const;

	const MaterialPtr& getMaterial() const;

	bool isCameraLens() const;

	bool isLight() const;

	void computeUVDifferential(const RayDifferential& ray);

	Fragment fragment;

	const Primitive* primitive;
};

typedef bool (*IntersectFilter)(const Primitive* p, const Ray& ray);

class Primitive {
public:
	virtual ~Primitive() = default;

	virtual bool intersectable() const {
		return true;
	}

	virtual void refine(PrimitiveList& refinedPrimitives) const {
		std::cerr << "unimplemented Primitive::refine" << std::endl;
		throw std::exception();
	}

	virtual bool intersect(const Ray& ray, float* epsilon,
		Intersection* intersection, IntersectFilter f = nullptr) const = 0;

	virtual bool occluded(const Ray& ray,
		IntersectFilter f = nullptr) const = 0;

	virtual BBox getAABB() const = 0;

	virtual const MaterialPtr& getMaterial() const {
		std::cerr << "unimplemented Primitive::getMaterial" << std::endl;
		throw std::exception();
	}

	virtual const AreaLight* getAreaLight() const {
		std::cerr << "unimplemented Primitive::getAreaLight" << std::endl;
		throw std::exception();
	}

	virtual bool isCameraLens() const {
		std::cerr << "unimplemented Primitive::isCameraLens" << std::endl;
		throw std::exception();
	}

	Primitive() {}
};

class InstancedPrimitive : public Primitive {
public:
	InstancedPrimitive(const Transform& toWorld,
		const Primitive* primitive);

	bool intersect(const Ray& ray, float* epsilon,
		Intersection* intersection, IntersectFilter f) const override;

	bool occluded(const Ray& ray, IntersectFilter f) const override;

	BBox getAABB() const override;

private:
	Transform mToWorld;
	const Primitive* mPrimitive;
};

class ParamSet;
class SceneCache;

Primitive* createInstance(const ParamSet& params, const SceneCache& sceneCache);

}

#endif // GOBLIN_PRIMITIVE_H