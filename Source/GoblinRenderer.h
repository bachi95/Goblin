#ifndef GOBLIN_RENDERER_H
#define GOBLIN_RENDERER_H

#include "GoblinScene.h"

namespace Goblin {
    class CameraSample;
    class Color;
    class Ray;
    class Sampler;
    class Renderer {
    public:
        Renderer(int maxRayDepth = 5);
        ~Renderer();

        void render(ScenePtr scene);
    private:
        Color Li(ScenePtr scene, const Ray& ray) const;
        Color specularReflect(const ScenePtr& scene, const Ray& ray, 
            float epsilon, const Intersection& intersection) const;
        Color specularRefract(const ScenePtr& scene, const Ray& ray, 
            float epsilon, const Intersection& intersection) const;

    private:
        CameraSample* mSamples;
        Sampler* mSampler;
        int mMaxRayDepth;
    };
}

#endif //GOBLIN_RENDERER_H