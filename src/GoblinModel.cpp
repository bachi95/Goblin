#include "GoblinBBox.h"
#include "GoblinBVH.h"
#include "GoblinModel.h"
#include "GoblinParamSet.h"
#include "GoblinScene.h"
#include "GoblinLight.h"

namespace Goblin {
std::vector<Model*> Model::refinedModels;

Model::Model(const Geometry* geometry, const MaterialPtr& material,
    const AreaLight* areaLight, bool isCameraLens):
    mGeometry(geometry), mMaterial(material), mAreaLight(areaLight),
    mIsCameraLens(isCameraLens) {
	if (!mGeometry->intersectable()) {
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
    GeometryList refinedGeometries;
    mGeometry->refine(refinedGeometries);
        
    Model* refined = new Model[refinedGeometries.size()];
    for (size_t i = 0; i < refinedGeometries.size(); ++i) {
        refined[i].init(refinedGeometries[i], mMaterial, getAreaLight());
    }
    refinedModels.push_back(refined);

    for (size_t i = 0; i < refinedGeometries.size(); ++i) {
        refinedPrimitives.push_back(&refined[i]);
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
    Primitive* model = new Model(geometry, material, areaLight,
        isCameraLens);
    Primitive::allocatedPrimitives.push_back(model);
	return model;
}

}
