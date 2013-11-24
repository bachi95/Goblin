#ifndef GOBLIN_MODEL_H
#define GOBLIN_MODEL_H
#include "GoblinPrimitive.h"

namespace Goblin {
    class Model : public Primitive {
    public:
        Model(const GeometryPtr& geometry, const MaterialPtr& material);
        bool intersectable() const;
        bool intersect(const Ray& ray);
        bool intersect(const Ray& ray, float* epsilon, 
            Intersection* intersection);
        BBox getAABB() const;
        void refine(PrimitiveList& refinedPrimitives);
        void collectRenderList(RenderList& rList, 
            const Matrix4& m = Matrix4::Identity);
        const MaterialPtr& getMaterial() const;
    private:
        GeometryPtr mGeometry;
        MaterialPtr mMaterial;
    };

    inline bool Model::intersectable() const { 
        return mGeometry->intersectable(); 
    }

    inline const MaterialPtr& Model::getMaterial() const {
        return mMaterial;
    }
}

#endif // GOBLIN_MODEL_H
