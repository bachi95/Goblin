#ifndef GOBLIN_MODEL_H
#define GOBLIN_MODEL_H
#include "GoblinPrimitive.h"

namespace Goblin {
    class Model : public Primitive {
    public:
        Model(const GeometryPtr& geometry);
        bool intersectable() const;
        bool intersect(const Ray& ray);
        void refine(PrimitiveList& refinedPrimitives);
        void collectRenderList(RenderList& rList, 
            const Matrix4& m = Matrix4::Identity);
    private:
        GeometryPtr mGeometry;
        //and Material.....
    };

    inline bool Model::intersectable() const { 
        return mGeometry->intersectable(); 
    }
}

#endif // GOBLIN_MODEL_H
