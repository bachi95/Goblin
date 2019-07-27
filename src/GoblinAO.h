#ifndef GOBLIN_AO_H
#define GOBLIN_AO_H

#include "GoblinRenderer.h"

namespace Goblin {
class AORenderer : public Renderer {
public:
    AORenderer(int samplePerPixel, int threadNum,
        int aoSample);
    ~AORenderer();
    Color Li(const ScenePtr& scene, const RayDifferential& ray,
        const Sample& sample, const RNG& rng,
        RenderingTLS* tls) const;
private:
    void querySampleQuota(const ScenePtr& scene,
        SampleQuota* sampleQuota);
private:
    int mAOSampleNum;
    SampleIndex mAOSampleIndex;
};

Renderer* createAO(const ParamSet& params);

}

#endif //GOBLIN_AO_H
