#include "GoblinBBox.h"
#include "GoblinBVH.h"
#include "GoblinModel.h"
#include "GoblinParamSet.h"
#include "GoblinScene.h"
#include "GoblinLight.h"

namespace Goblin {

Model::Model(const Geometry* geometry, const MaterialPtr& material,
    const AreaLight* areaLight, bool isCameraLens):
    mGeometry(geometry), mMaterial(material), mAreaLight(areaLight),
    mIsCameraLens(isCameraLens) {
	if (!mGeometry->intersectable()) {
		GeometryList refinedGeometries;
		mGeometry->refine(refinedGeometries);
		mRefinedModels.reserve(refinedGeometries.size());
		for (size_t i = 0; i < refinedGeometries.size(); ++i) {
			mRefinedModels.emplace_back(
				refinedGeometries[i], mMaterial, getAreaLight(), mIsCameraLens);
		}
		PrimitiveList primitives;
		primitives.push_back(this);
		mBVH.reset(new BVH(primitives, 1, "equal_count"));
	}
}

bool Model::occluded(const Ray& ray, IntersectFilter f) const {
	if (mBVH) {
		return mBVH->occluded(ray, f);
	} else {
		if (f != nullptr && !f(this, ray)) {
			return false;
		}
		return mGeometry->occluded(ray);
	}
}

bool Model::intersect(const Ray& ray, float* epsilon, 
    Intersection* intersection, IntersectFilter f) const {
	if (mBVH) {
		return mBVH->intersect(ray, epsilon, intersection, f);
	} else {
		if (f != nullptr && !f(this, ray)) {
			return false;
		}
		bool hit = mGeometry->intersect(ray, epsilon,
			&intersection->fragment);
		if (hit) {
			intersection->primitive = this;
		}
		return hit;
	}
}

BBox Model::getAABB() const {
    return mGeometry->getObjectBound();
}

void Model::refine(PrimitiveList& refinedPrimitives) const {
    for (size_t i = 0; i < mRefinedModels.size(); ++i) {
		refinedPrimitives.push_back(&mRefinedModels[i]);
    }
}

Primitive* createModel(const ParamSet& params,
    const SceneCache& sceneCache) {

	std::string geoName = params.getString("geometry");
    const Geometry* geometry = sceneCache.getGeometry(geoName);

	std::string materialName = params.getString("material");
    MaterialPtr material = sceneCache.getMaterial(materialName);

    const AreaLight* areaLight = nullptr;
    if (params.hasString("area_light")) {
        areaLight = sceneCache.getAreaLight(params.getString("area_light"));
    }
    bool isCameraLens = params.getBool("is_camera_lens");
	return new Model(geometry, material, areaLight,
		isCameraLens);
}

}
