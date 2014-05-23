#ifndef GOBLIN_PATHTRACER_H
#define GOBLIN_PATHTRACER_H

#include "GoblinRenderer.h"

namespace Goblin {
    class PathTracer : public Renderer {
    public:
        PathTracer(const RenderSetting& setting);
        ~PathTracer();
        Color Li(const ScenePtr& scene, const Ray& ray,
            const Sample& sample, const RNG& rng) const;
    private:
        void querySampleQuota(const ScenePtr& scene,
            SampleQuota* sampleQuota);
    private:
        BSDFSampleIndex* mPathSampleIndexes;
    };
}

#endif //GOBLIN_PATHTRACER_H
