#ifndef GOBLIN_RENDERER_H
#define GOBLIN_RENDERER_H

#include "GoblinMaterial.h"
#include "GoblinScene.h"
#include "GoblinSampler.h"
#include "GoblinThreadPool.h"

namespace Goblin {
    class Color;
    class ImageTile;
    class Ray;
    class ParamSet;
    class Renderer;
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

    class RenderTask : public Task {
    public:
        RenderTask(ImageTile* tile, Renderer* mRenderer,
            const CameraPtr& camera, const ScenePtr& scene,
            const SampleQuota& sampleQuota, 
            int samplePerPixel);
        ~RenderTask();
        void run();
    private:
        ImageTile* mTile;
        Renderer* mRenderer;
        const CameraPtr& mCamera;
        const ScenePtr& mScene;
        const SampleQuota& mSampleQuota;
        int mSamplePerPixel;
        RNG* mRNG;
    };

    class Renderer {
    public:
        Renderer(const RenderSetting& setting);
        virtual ~Renderer();
        void render(const ScenePtr& scene);
        virtual Color Li(const ScenePtr& scene, const Ray& ray, 
            const Sample& sample, const RNG& rng) const = 0;

    protected:
        Color singleSampleLd(const ScenePtr& scene, const Ray& ray,
            float epsilon, const Intersection& intersection, 
            const Sample& sample, 
            const LightSample& lightSample,
            const BSDFSample& bsdfSample,
            float pickLightSample,
            BSDFType type = BSDFAll) const;

        Color multiSampleLd(const ScenePtr& scene, const Ray& ray,
            float epsilon, const Intersection& intersection, 
            const Sample& sample, const RNG& rng,
            LightSampleIndex* lightSampleIndexes = NULL,
            BSDFSampleIndex* bsdfSampleIndexes = NULL,
            BSDFType type = BSDFAll) const;

        Color estimateLd(const ScenePtr& scene, const Vector3& wo,
            float epsilon, const Intersection& intersection, 
            const Light* light, const LightSample& ls,
            const BSDFSample& bs, BSDFType type = BSDFAll) const;

        Color specularReflect(const ScenePtr& scene, const Ray& ray, 
            float epsilon, const Intersection& intersection,
            const Sample& sample, const RNG& rng) const;

        Color specularRefract(const ScenePtr& scene, const Ray& ray, 
            float epsilon, const Intersection& intersection,
            const Sample& sample, const RNG& rng) const;
    private:
        virtual void querySampleQuota(const ScenePtr& scene, 
            SampleQuota* sampleQuota) = 0;

    protected:
        LightSampleIndex* mLightSampleIndexes;
        BSDFSampleIndex* mBSDFSampleIndexes;
        SampleIndex* mPickLightSampleIndexes;
        CDF1D* mPowerDistribution;
        RenderSetting mSetting;
    };
}

#endif //GOBLIN_RENDERER_H
