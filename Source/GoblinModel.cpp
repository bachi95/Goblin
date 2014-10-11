#include "GoblinBBox.h"
#include "GoblinBVH.h"
#include "GoblinModel.h"
#include "GoblinParamSet.h"
#include "GoblinScene.h"

namespace Goblin {
    Model::Model(const GeometryPtr& geometry, const MaterialPtr& material,
        const AreaLight* areaLight): 
        mGeometry(geometry), mMaterial(material), mAreaLight(areaLight) {}

    bool Model::intersect(const Ray& ray) {
        return mGeometry->intersect(ray);
    }

    bool Model::intersect(const Ray& ray, float* epsilon, 
        Intersection* intersection) {
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

    void Model::refine(PrimitiveList& refinedPrimitives) {
        GeometryList refinedGeometries;
        mGeometry->refine(refinedGeometries);
        for(size_t i = 0; i < refinedGeometries.size(); ++i) {
            PrimitivePtr primitive(new Model(refinedGeometries[i], mMaterial,
                getAreaLight()));
            refinedPrimitives.push_back(primitive);
        }
    }

    void Model::collectRenderList(RenderList& rList, 
        const Matrix4& m) {
        rList.push_back(Renderable(m, mGeometry));
    }


    Primitive* ModelPrimitiveCreator::create(const ParamSet& params,
        const SceneCache& sceneCache) const {

        string geoName = params.getString("geometry");
        GeometryPtr geometry = sceneCache.getGeometry(geoName);

        string materialName = params.getString("material");
        MaterialPtr material = sceneCache.getMaterial(materialName);

        Primitive* model = new Model(geometry, material);

        if(!model->intersectable()) {
            std::vector<PrimitivePtr> primitives;
            primitives.push_back(PrimitivePtr(model));
            Primitive* aggregate = new BVH(primitives, 1, "equal_count");
            return aggregate;
        } else {
            return model;
        }
    }

}
