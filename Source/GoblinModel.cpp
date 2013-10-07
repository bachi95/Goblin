#include "GoblinModel.h"

namespace Goblin {
    Model::Model(const GeometryPtr& geometry) : mGeometry(geometry) {}

    bool Model::intersect(const Ray& ray) {
        return mGeometry->intersect(ray);
    }

    void Model::collectRenderList(RenderList& rList, 
        const Matrix4& m) {
        rList.push_back(Renderable(m, mGeometry));
    }
}