#ifndef GOBLIN_PATHTRACER_H
#define GOBLIN_PATHTRACER_H

#include "GoblinRenderer.h"

namespace Goblin {
    class PathTracer : public Renderer {
    public:
        PathTracer(const ParamSet& setting);
        ~PathTracer();
        Color Li(const ScenePtr& scene, const RayDifferential& ray,
            const Sample& sample, const RNG& rng,
            WorldDebugData* debugData) const;
    private:
        // evaluate index-matched material attenuation along the ray
        Color evalAttenuation(const ScenePtr& scene, const Ray& ray,
            const BSDFSample& bs) const;
        
        void querySampleQuota(const ScenePtr& scene,
            SampleQuota* sampleQuota);
    private:
        int mMaxRayDepth;
    };
}

#endif //GOBLIN_PATHTRACER_H
