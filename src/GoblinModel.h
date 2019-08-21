#ifndef GOBLIN_MODEL_H
#define GOBLIN_MODEL_H
#include "GoblinPrimitive.h"
#include "GoblinLight.h"

namespace Goblin {

class BVH;

class Model : public Primitive {
public:
	Model(const Geometry* geometry, const MaterialPtr& material,
		const AreaLight* areaLight, bool isCameraLens);

	bool intersectable() const override {
		return mBVH ? true : mGeometry->intersectable();
	}

	bool intersect(const Ray& ray, float* epsilon,
		Intersection* intersection, IntersectFilter f) const override;

	bool occluded(const Ray& ray, IntersectFilter f) const override;

	bool isCameraLens() const override {
		return mIsCameraLens;
	}

    BBox getAABB() const override;

	const MaterialPtr& getMaterial() const override {
		return mMaterial;
	}

	const AreaLight* getAreaLight() const override {
		return mAreaLight;
	}

    void refine(PrimitiveList& refinedPrimitives) const override;

private:
    const Geometry* mGeometry;
    MaterialPtr mMaterial;
    const AreaLight* mAreaLight;
    bool mIsCameraLens;
    std::vector<Model> mRefinedModels;
	std::unique_ptr<BVH> mBVH;
};

class ParamSet;
class SceneCache;

Primitive* createModel(const ParamSet& params,
	const SceneCache& sceneCache);

}

#endif // GOBLIN_MODEL_H