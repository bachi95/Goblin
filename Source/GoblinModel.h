#ifndef GOBLIN_MODEL_H
#define GOBLIN_MODEL_H
#include "GoblinFactory.h"
#include "GoblinPrimitive.h"
#include "GoblinLight.h"

namespace Goblin {
    class Model : public Primitive {
    public:
        bool intersectable() const;
        bool intersect(const Ray& ray, IntersectFilter f) const;
        bool intersect(const Ray& ray, float* epsilon, 
            Intersection* intersection, IntersectFilter f) const;

        BBox getAABB() const;
        const MaterialPtr& getMaterial() const;
        const AreaLight* getAreaLight() const;

        void refine(PrimitiveList& refinedPrimitives) const;
        void collectRenderList(RenderList& rList, 
            const Matrix4& m = Matrix4::Identity) const ;
        static void clearRefinedModels();
    private:
        Model(const Geometry* geometry, const MaterialPtr& material,
            const AreaLight* areaLight = NULL);
        // For chunk allocation refined models, 
        // allocate first then set geometry
        Model() {};
        void init(const Geometry* geometry, const MaterialPtr& material,
            const AreaLight* areaLight);
    private:
        const Geometry* mGeometry;
        MaterialPtr mMaterial;
        const AreaLight* mAreaLight;
        // used to keep the refined models generated by refine method
        static vector<Model*> refinedModels;
    friend class ModelPrimitiveCreator;
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

    inline void Model::init(const Geometry* geometry, 
        const MaterialPtr& material, const AreaLight* areaLight) {
        mGeometry = geometry;
        mMaterial = material;
        mAreaLight = areaLight;
    }

    inline void Model::clearRefinedModels() {
        for(size_t i = 0; i < refinedModels.size(); ++i) {
            delete [] refinedModels[i];
            refinedModels[i] = NULL;
        }
    }
    
    class ParamSet;
    class SceneCache;

    class ModelPrimitiveCreator : public 
        Creator<Primitive , const ParamSet&, const SceneCache&> {
    public:
        Primitive* create(const ParamSet& params, 
            const SceneCache& sceneCache) const;
    };

}

#endif // GOBLIN_MODEL_H
