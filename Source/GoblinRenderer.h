#ifndef GOBLIN_RENDERER_H
#define GOBLIN_RENDERER_H

#include "GoblinMaterial.h"
#include "GoblinScene.h"
#include "GoblinSampler.h"

namespace Goblin {
    class Color;
    class Ray;
    class ParamSet;
    struct BSDFSample;
    struct LightSample;
    struct BSDFSampleIndex;
    struct LightSampleIndex;

    struct RenderSetting {
        RenderSetting(): 
            samplePerPixel(1), maxRayDepth(5) {}
        int samplePerPixel;
        int maxRayDepth;
    };

    class Renderer {
    public:
        Renderer(const RenderSetting& setting);
        virtual ~Renderer();
        void render(const ScenePtr& scene);

    protected:
        virtual Color Li(const ScenePtr& scene, const Ray& ray, 
            const Sample& sample) const = 0;

        Color singleSampleLd(const ScenePtr& scene, const Ray& ray,
            float epsilon, const Intersection& intersection, 
            const Sample& sample, 
            const LightSample& lightSample,
            const BSDFSample& bsdfSample,
            float pickLightSample,
            BSDFType type = BSDFAll) const;

        Color multiSampleLd(const ScenePtr& scene, const Ray& ray,
            float epsilon, const Intersection& intersection, 
            const Sample& sample, LightSampleIndex* lightSampleIndexes = NULL,
            BSDFSampleIndex* bsdfSampleIndexes = NULL,
            BSDFType type = BSDFAll) const;

        Color estimateLd(const ScenePtr& scene, const Vector3& wo,
            float epsilon, const Intersection& intersection, 
            const Light* light, const LightSample& ls,
            const BSDFSample& bs, BSDFType type = BSDFAll) const;

        Color specularReflect(const ScenePtr& scene, const Ray& ray, 
            float epsilon, const Intersection& intersection,
            const Sample& sample) const;

        Color specularRefract(const ScenePtr& scene, const Ray& ray, 
            float epsilon, const Intersection& intersection,
            const Sample& sample) const;
    private:
        virtual void querySampleQuota(const ScenePtr& scene, 
            Sampler* sampler) = 0;

    protected:
        LightSampleIndex* mLightSampleIndexes;
        BSDFSampleIndex* mBSDFSampleIndexes;
        SampleIndex* mPickLightSampleIndexes;
        Sample* mSamples;
        Sampler* mSampler;
        CDF1D* mPowerDistribution;
        RenderSetting mSetting;
    };
}

#endif //GOBLIN_RENDERER_H
