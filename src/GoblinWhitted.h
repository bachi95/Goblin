#ifndef GOBLIN_WHITTED_H
#define GOBLIN_WHITTED_H

#include "GoblinRenderer.h"

namespace Goblin {

class WhittedRenderer : public Renderer {
public:
    WhittedRenderer(int samplePerPixel, int threadNum, 
        int maxRayDepth, int bssrdfSampleNum);

    ~WhittedRenderer() = default;

    Color Li(const ScenePtr& scene, const RayDifferential& ray, 
        const Sample& sample, const RNG& rng,
        RenderingTLS* tls) const override;

private:
    void querySampleQuota(const ScenePtr& scene, 
        SampleQuota* sampleQuota) override;

private:
    int mMaxRayDepth;
    int mBssrdfSampleNum;
};

Renderer* createWhitted(const ParamSet& params);

}

#endif // GOBLIN_WHITTED_H