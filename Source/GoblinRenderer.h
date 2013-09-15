#ifndef GOBLIN_RENDERER_H
#define GOBLIN_RENDERER_H

namespace Goblin {
    class CameraSample;
    class Color;
    class Scene;
    class Ray;
    class Sampler;
    class Renderer {
    public:
        Renderer();
        ~Renderer();

        void render(Scene* scene);
        Color Li(Scene* scene, const Ray& ray);
    private:
        CameraSample* mSamples;
        Sampler* mSampler;
    };
}

#endif //GOBLIN_RENDERER_H