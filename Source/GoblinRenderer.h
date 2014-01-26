#ifndef GOBLIN_RENDERER_H
#define GOBLIN_RENDERER_H

#include "GoblinScene.h"
#include "GoblinSampler.h"

namespace Goblin {
    class Sample;
    class Color;
    class Ray;
    class Sampler;
    class ParamSet;

    struct RenderSetting {
        RenderSetting(): 
            samplePerPixel(1), maxRayDepth(5) {}
        int samplePerPixel;
        int maxRayDepth;
    };

    class Renderer {
    public:
        Renderer(const RenderSetting& setting);
        ~Renderer();
        void render(const ScenePtr& scene);
    private:
        Color Li(const ScenePtr& scene, const Ray& ray) const;
        Color specularReflect(const ScenePtr& scene, const Ray& ray, 
            float epsilon, const Intersection& intersection) const;
        Color specularRefract(const ScenePtr& scene, const Ray& ray, 
            float epsilon, const Intersection& intersection) const;

        void querySampleQuota(const ScenePtr& scene, Sampler* sampler);
    private:
        Sample* mSamples;
        Sampler* mSampler;
        int mMaxRayDepth;
        RenderSetting mSetting;
    };
}

#endif //GOBLIN_RENDERER_H