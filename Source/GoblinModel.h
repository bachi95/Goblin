#ifndef GOBLIN_MODEL_H
#define GOBLIN_MODEL_H
#include "GoblinPrimitive.h"
#include "GoblinLight.h"

namespace Goblin {
    class Model : public Primitive {
    public:
        Model(const GeometryPtr& geometry, const MaterialPtr& material,
            const AreaLight* areaLight = NULL);
        bool intersectable() const;
        bool intersect(const Ray& ray);
        bool intersect(const Ray& ray, float* epsilon, 
            Intersection* intersection);
        BBox getAABB() const;
        const MaterialPtr& getMaterial() const;
        const AreaLight* getAreaLight() const;

        void refine(PrimitiveList& refinedPrimitives);
        void collectRenderList(RenderList& rList, 
            const Matrix4& m = Matrix4::Identity);
    private:
        GeometryPtr mGeometry;
        MaterialPtr mMaterial;
        const AreaLight* mAreaLight;
    };

    inline bool Model::intersectable() const { 
        return mGeometry->intersectable(); 
    }

    inline const MaterialPtr& Model::getMaterial() const {
        return mMaterial;
    }

    inline const AreaLight* Model::getAreaLight() const {
        return mAreaLight;
    }
}

#endif // GOBLIN_MODEL_H
