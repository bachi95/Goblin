#ifndef GOBLIN_AO_H
#define GOBLIN_AO_H

#include "GoblinFactory.h"
#include "GoblinRenderer.h"

namespace Goblin {
    class AORenderer : public Renderer {
    public:
        AORenderer(int samplePerPixel = 1, int threadNum = 1, 
            int aoSample = 25);
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

    class AORendererCreator : public 
        Creator<Renderer , const ParamSet&> {
    public:
        Renderer* create(const ParamSet& params) const;
    };
}

#endif //GOBLIN_AO_H
