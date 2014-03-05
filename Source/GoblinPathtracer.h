#ifndef GOBLIN_PATHTRACER_H
#define GOBLIN_PATHTRACER_H

#include "GoblinRenderer.h"

namespace Goblin {
    class PathTracer : public Renderer {
    public:
        PathTracer(const RenderSetting& setting);
        ~PathTracer();
        Color Li(const ScenePtr& scene, const Ray& ray,
            const Sample& sample) const;
    private:
        void querySampleQuota(const ScenePtr& scene,
            Sampler* sampler);
    private:
        BSDFSampleIndex* mPathSampleIndexes;
    };
}

#endif //GOBLIN_PATHTRACER_H
