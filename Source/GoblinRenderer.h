#ifndef GOBLIN_RENDERER_H
#define GOBLIN_RENDERER_H

#include "GoblinMaterial.h"
#include "GoblinRay.h"
#include "GoblinScene.h"
#include "GoblinSampler.h"
#include "GoblinThreadPool.h"

namespace Goblin {
    class Color;
    class ImageTile;
    class ParamSet;
    class Renderer;
    struct BSDFSample;
    struct LightSample;
    struct BSDFSampleIndex;
    struct LightSampleIndex;

    class RenderProgress {
    public:
        RenderProgress(int taskNum);
        void reset();
        void update();
    private:
        boost::mutex mUpdateMutex;
        int mFinishedNum, mTasksNum;
    };

    class RenderTask : public Task {
    public:
        RenderTask(Renderer* mRenderer, const CameraPtr& camera,
            const ScenePtr& scene, const SampleRange& sampleRange,
            const SampleQuota& sampleQuota, int samplePerPixel, 
            RenderProgress* renderProgress);
        ~RenderTask();
        void run(TLSPtr& tls);

    protected:
        Renderer* mRenderer;
        const CameraPtr& mCamera;
        const ScenePtr& mScene;
        const SampleRange& mSampleRange;
        const SampleQuota& mSampleQuota;
        int mSamplePerPixel;
        RenderProgress* mRenderProgress;
        RNG* mRNG;
    };

    class Renderer {
    public:
        Renderer(int samplePerPixel = 1, int threadNum = 1);
        virtual ~Renderer();
        virtual void render(const ScenePtr& scene);
        virtual Color Li(const ScenePtr& scene, const RayDifferential& ray, 
            const Sample& sample, const RNG& rng,
            RenderingTLS* tls = NULL) const = 0;
        // volume in scatter and emission contribution
        Color Lv(const ScenePtr& scene, const Ray& ray, const RNG& rng) const;
        Color Lsubsurface(const ScenePtr& scene,
            const Intersection& intersection, const Vector3& wo,
            const Sample& sample, const BSSRDFSampleIndex* bssrdfSampleIndex,
            RenderingTLS* tls = NULL) const;
        Color transmittance(const ScenePtr& scene, const Ray& ray) const;

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

        Color singleSampleIrradiance(const ScenePtr& scene,
            float epsilon, const Intersection& intersection, 
            const LightSample& lightSample, float pickLightSample) const;

        Color estimateIrradiance(const ScenePtr& scene,
            float epsilon, const Intersection& intersection, 
            const Light* light, const LightSample& ls) const;

        Color specularReflect(const ScenePtr& scene, 
            const RayDifferential& ray, 
            float epsilon, const Intersection& intersection,
            const Sample& sample, const RNG& rng) const;

        Color specularRefract(const ScenePtr& scene, 
            const RayDifferential& ray, 
            float epsilon, const Intersection& intersection,
            const Sample& sample, const RNG& rng) const;

        void getSampleRanges(const Film* film,
            vector<SampleRange>& sampleRanges) const;

        void drawDebugData(const DebugData& debugData,
            const CameraPtr& camera) const;

    private:
        virtual void querySampleQuota(const ScenePtr& scene, 
            SampleQuota* sampleQuota) = 0;

        Color LbssrdfSingle(const ScenePtr& scene, const Fragment& fragment, 
            const BSSRDF* bssrdf, const Vector3& wo, const Sample& sample, 
            const BSSRDFSampleIndex* bssrdfSampleIndex,
            RenderingTLS* tls = NULL) const;

        Color LbssrdfDiffusion(const ScenePtr& scene, const Fragment& fragment, 
            const BSSRDF* bssrdf, const Vector3& wo, const Sample& sample, 
            const BSSRDFSampleIndex* bssrdfSampleIndex,
            RenderingTLS* tls = NULL ) const;


    protected:
        LightSampleIndex* mLightSampleIndexes;
        BSDFSampleIndex* mBSDFSampleIndexes;
        SampleIndex* mPickLightSampleIndexes;
        BSSRDFSampleIndex mBSSRDFSampleIndex;
        CDF1D* mPowerDistribution;
        int mSamplePerPixel;
        int mThreadNum;
    };
}

#endif //GOBLIN_RENDERER_H
