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

    enum RenderMethod {
        Whitted,
        PathTracing
    };

    struct RenderSetting {
        RenderSetting(): samplePerPixel(1), threadNum(1), maxRayDepth(5), 
            bssrdfSampleNum(4), method(PathTracing) {}
        int samplePerPixel;
        int threadNum;
        int maxRayDepth;
        int bssrdfSampleNum;
        RenderMethod method;
    };

    class WorldDebugData {
    public:
        WorldDebugData() {}
        void addRay(const Ray& ray, Color c = Color::White);
        void addPoint(const Vector3& point, Color c = Color::White);
        const vector<pair<Ray, Color> >& getRays() const;
        const vector<pair<Vector3, Color> >& getPoints() const;
    private:
        vector<pair<Ray, Color> > mRays;
        vector<pair<Vector3, Color> > mPoints;
    };

    inline void WorldDebugData::addRay(const Ray& line, Color c) { 
        mRays.push_back(pair<Ray, Color>(line, c)); 
    }

    inline void WorldDebugData::addPoint(const Vector3& point, Color c) { 
        mPoints.push_back(pair<Vector3, Color>(point, c)); 
    }

    inline const vector<pair<Ray, Color> >& WorldDebugData::getRays() const {
        return mRays;
    }

    inline const vector<pair<Vector3, Color> >& WorldDebugData::getPoints() 
        const {
        return mPoints;
    }


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
        RenderTask(ImageTile* tile, Renderer* mRenderer,
            const CameraPtr& camera, const ScenePtr& scene,
            const SampleQuota& sampleQuota, int samplePerPixel, 
            RenderProgress* renderProgress);
        ~RenderTask();
        void run();
    private:
        ImageTile* mTile;
        Renderer* mRenderer;
        const CameraPtr& mCamera;
        const ScenePtr& mScene;
        const SampleQuota& mSampleQuota;
        int mSamplePerPixel;
        RenderProgress* mRenderProgress;
        RNG* mRNG;
    };

    class Renderer {
    public:
        Renderer(const RenderSetting& setting);
        virtual ~Renderer();
        void render(const ScenePtr& scene);
        virtual Color Li(const ScenePtr& scene, const Ray& ray, 
            const Sample& sample, const RNG& rng,
            WorldDebugData* debugData = NULL) const = 0;
        // volume in scatter and emission contribution
        Color Lv(const ScenePtr& scene, const Ray& ray, const RNG& rng) const;
        Color Lsubsurface(const ScenePtr& scene,
            const Intersection& intersection, const Vector3& wo,
            const Sample& sample, const BSSRDFSampleIndex* bssrdfSampleIndex,
            WorldDebugData* debugData = NULL) const;
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

        Color specularReflect(const ScenePtr& scene, const Ray& ray, 
            float epsilon, const Intersection& intersection,
            const Sample& sample, const RNG& rng) const;

        Color specularRefract(const ScenePtr& scene, const Ray& ray, 
            float epsilon, const Intersection& intersection,
            const Sample& sample, const RNG& rng) const;

    private:
        virtual void querySampleQuota(const ScenePtr& scene, 
            SampleQuota* sampleQuota) = 0;

        Color LbssrdfSingle(const ScenePtr& scene, const Fragment& fragment, 
            const BSSRDF* bssrdf, const Vector3& wo, const Sample& sample, 
            const BSSRDFSampleIndex* bssrdfSampleIndex) const;

        Color LbssrdfDiffusion(const ScenePtr& scene, const Fragment& fragment, 
            const BSSRDF* bssrdf, const Vector3& wo, const Sample& sample, 
            const BSSRDFSampleIndex* bssrdfSampleIndex) const;


    protected:
        LightSampleIndex* mLightSampleIndexes;
        BSDFSampleIndex* mBSDFSampleIndexes;
        SampleIndex* mPickLightSampleIndexes;
        BSSRDFSampleIndex mBSSRDFSampleIndex;
        CDF1D* mPowerDistribution;
        RenderSetting mSetting;
    };
}

#endif //GOBLIN_RENDERER_H
