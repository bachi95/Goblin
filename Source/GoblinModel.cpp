#include "GoblinBBox.h"
#include "GoblinBVH.h"
#include "GoblinModel.h"
#include "GoblinParamSet.h"
#include "GoblinScene.h"
#include "GoblinLight.h"

namespace Goblin {
    vector<Model*> Model::refinedModels;

    Model::Model(const Geometry* geometry, const MaterialPtr& material,
        const AreaLight* areaLight): 
        mGeometry(geometry), mMaterial(material), mAreaLight(areaLight) {}

    bool Model::intersect(const Ray& ray, IntersectFilter f) const {
        if(f != NULL && !f(this, ray)) {
            return false;
        }
        return mGeometry->intersect(ray);
    }

    bool Model::intersect(const Ray& ray, float* epsilon, 
        Intersection* intersection, IntersectFilter f) const {
        if(f != NULL && !f(this, ray)) {
            return false;
        }
        bool hit = mGeometry->intersect(ray, epsilon, 
            &intersection->fragment);
        if(hit) {
            intersection->primitive = this;
        }
        return hit;
    }

    BBox Model::getAABB() const {
        return mGeometry->getObjectBound();
    }

    void Model::refine(PrimitiveList& refinedPrimitives) const {
        GeometryList refinedGeometries;
        mGeometry->refine(refinedGeometries);
        
        Model* refined = new Model[refinedGeometries.size()];
        for(size_t i = 0; i < refinedGeometries.size(); ++i) {
            refined[i].init(refinedGeometries[i], mMaterial, getAreaLight());
        }
        refinedModels.push_back(refined);

        for(size_t i = 0; i < refinedGeometries.size(); ++i) {
            refinedPrimitives.push_back(&refined[i]);
        }
    }

    void Model::collectRenderList(RenderList& rList, 
        const Matrix4& m) const {
        rList.push_back(Renderable(m, mGeometry));
    }


    Primitive* ModelPrimitiveCreator::create(const ParamSet& params,
        const SceneCache& sceneCache) const {

        string geoName = params.getString("geometry");
        const Geometry* geometry = sceneCache.getGeometry(geoName);

        string materialName = params.getString("material");
        MaterialPtr material = sceneCache.getMaterial(materialName);

        const AreaLight* areaLight = NULL;
        if(params.hasString("area_light")) {
            areaLight = sceneCache.getAreaLight(params.getString("area_light"));
        }
        Primitive* model = new Model(geometry, material, areaLight);
        Primitive::allocatedPrimitives.push_back(model);

        if(!model->intersectable()) {
            PrimitiveList primitives;
            primitives.push_back(model);
            Primitive* aggregate = new BVH(primitives, 1, "equal_count");
            Primitive::allocatedPrimitives.push_back(aggregate);
            return aggregate;
        } else {
            return model;
        }
    }

}
