#ifndef GOBLIN_WHITTED_H
#define GOBLIN_WHITTED_H

#include "GoblinRenderer.h"

namespace Goblin {
class WhittedRenderer : public Renderer {
public:
    WhittedRenderer(int samplePerPixel = 1, int threadNum = 1, 
        int maxRayDepth = 5, int bssrdfSampleNum = 4);
    ~WhittedRenderer();
    Color Li(const ScenePtr& scene, const RayDifferential& ray, 
        const Sample& sample, const RNG& rng,
        RenderingTLS* tls) const;
private:
    void querySampleQuota(const ScenePtr& scene, 
        SampleQuota* sampleQuota);
private:
    int mMaxRayDepth;
    int mBssrdfSampleNum;
};

Renderer* createWhitted(const ParamSet& params);

}

#endif //GOBLIN_WHITTED_H
