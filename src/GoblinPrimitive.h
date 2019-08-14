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
typedef std::shared_ptr<Primitive> PrimitivePtr;
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
	virtual ~Primitive() {}

	virtual bool intersectable() const;

	virtual void refine(PrimitiveList& refinedPrimitives) const;

	virtual bool intersect(const Ray& ray, float* epsilon,
		Intersection* intersection, IntersectFilter f = nullptr) const = 0;

	virtual bool occluded(const Ray& ray,
		IntersectFilter f = nullptr) const = 0;

	virtual BBox getAABB() const = 0;

	virtual const MaterialPtr& getMaterial() const;

	virtual const AreaLight* getAreaLight() const;

	virtual bool isCameraLens() const;

	static void clearAllocatedPrimitives();

	static std::vector<Primitive*> allocatedPrimitives;

	Primitive() {}
};

inline bool Primitive::intersectable() const { return true; }

inline void Primitive::refine(PrimitiveList& refinedPrimitives) const {
	std::cerr << "unimplemented Primitive::refine" << std::endl;
	throw std::exception();
}

// TODO meh......instance should have the right to do 
//material override, add it in later...
inline const MaterialPtr& Primitive::getMaterial() const {
	std::cerr << "unimplemented Primitive::getMaterial" << std::endl;
	throw std::exception();
}

inline const AreaLight* Primitive::getAreaLight() const {
	std::cerr << "unimplemented Primitive::getAreaLight" << std::endl;
	throw std::exception();
}

inline bool Primitive::isCameraLens() const {
	std::cerr << "unimplemented Primitive::isCameraLens" << std::endl;
	throw std::exception();
}

inline void Primitive::clearAllocatedPrimitives() {
	for (size_t i = 0; i < allocatedPrimitives.size(); ++i) {
		delete allocatedPrimitives[i];
		allocatedPrimitives[i] = nullptr;
	}
}

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