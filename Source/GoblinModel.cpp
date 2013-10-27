#include "GoblinModel.h"

namespace Goblin {
    Model::Model(const GeometryPtr& geometry) : mGeometry(geometry) {}

    bool Model::intersect(const Ray& ray) {
        return mGeometry->intersect(ray);
    }

    bool Model::intersect(const Ray& ray, float* epsilon, 
        Intersection* intersection) {
        return mGeometry->intersect(ray, epsilon, intersection);
    }

    void Model::refine(PrimitiveList& refinedPrimitives) {
        GeometryList refinedGeometries;
        mGeometry->refine(refinedGeometries);
        for(size_t i = 0; i < refinedGeometries.size(); ++i) {
            PrimitivePtr primitive(new Model(refinedGeometries[i]));
            refinedPrimitives.push_back(primitive);
        }
    }

    void Model::collectRenderList(RenderList& rList, 
        const Matrix4& m) {
        rList.push_back(Renderable(m, mGeometry));
    }
}
