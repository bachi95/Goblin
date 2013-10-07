#ifndef GOBLIN_MODEL_H
#define GOBLIN_MODEL_H
#include "GoblinPrimitive.h"

namespace Goblin {
    class Model : public Primitive {
    public:
        Model(const GeometryPtr& geometry);
        bool intersect(const Ray& ray);
        void collectRenderList(RenderList& rList, 
            const Matrix4& m = Matrix4::Identity);
    private:
        GeometryPtr mGeometry;
        //and Material.....
    };
}

#endif // GOBLIN_MODEL_H