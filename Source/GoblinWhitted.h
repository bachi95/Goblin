#ifndef GOBLIN_WHITTED_H
#define GOBLIN_WHITTED_H

#include "GoblinRenderer.h"

namespace Goblin {
    class WhittedRenderer : public Renderer {
    public:
        WhittedRenderer(const ParamSet& setting);
        ~WhittedRenderer();
        Color Li(const ScenePtr& scene, const Ray& ray, 
            const Sample& sample, const RNG& rng,
            WorldDebugData* debugData) const;
    private:
        void querySampleQuota(const ScenePtr& scene, 
            SampleQuota* sampleQuota);
    private:
        int mMaxRayDepth;
    };
}

#endif //GOBLIN_WHITTED_H
