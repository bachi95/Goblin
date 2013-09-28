#include "GoblinModel.h"

namespace Goblin {
    Model::Model(const GeometryPtr& geometry) : mGeometry(geometry) {}

    void Model::collectRenderList(RenderList& rList, 
        const Matrix4& m) {
        rList.push_back(Renderable(m, mGeometry));
    }
}