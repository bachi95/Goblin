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
        Renderer();
        ~Renderer();

        void render(ScenePtr scene);
        Color Li(ScenePtr scene, const Ray& ray);
    private:
        CameraSample* mSamples;
        Sampler* mSampler;
    };
}

#endif //GOBLIN_RENDERER_H