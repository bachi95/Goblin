#ifndef GOBLIN_PATHTRACER_H
#define GOBLIN_PATHTRACER_H

#include "GoblinFactory.h"
#include "GoblinRenderer.h"

namespace Goblin {
class PathTracer : public Renderer {
public:
    PathTracer(int samplePerPixel = 1, int threadNum = 1, 
        int maxRayDepth = 5, int bssrdfSampleNum = 4);
    ~PathTracer();
    Color Li(const ScenePtr& scene, const RayDifferential& ray,
        const Sample& sample, const RNG& rng,
        RenderingTLS* tls) const;
private:
    // evaluate index-matched material attenuation along the ray
    Color evalAttenuation(const ScenePtr& scene, const Ray& ray,
        const BSDFSample& bs) const;
        
    void querySampleQuota(const ScenePtr& scene,
        SampleQuota* sampleQuota);
private:
    int mMaxRayDepth;
    int mBssrdfSampleNum;
};

class PathTracerCreator : public 
    Creator<Renderer , const ParamSet&> {
public:
    Renderer* create(const ParamSet& params) const;
};
}

#endif //GOBLIN_PATHTRACER_H
