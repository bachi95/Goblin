#include "GoblinModel.h"
#include "GoblinBBox.h"

namespace Goblin {
    Model::Model(const GeometryPtr& geometry, const MaterialPtr& material): 
        mGeometry(geometry), mMaterial(material) {}

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
            PrimitivePtr primitive(new Model(refinedGeometries[i], mMaterial));
            refinedPrimitives.push_back(primitive);
        }
    }

    void Model::collectRenderList(RenderList& rList, 
        const Matrix4& m) {
        rList.push_back(Renderable(m, mGeometry));
    }
}
