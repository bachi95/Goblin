#ifndef GOBLIN_RENDERER_H
#define GOBLIN_RENDERER_H

namespace Goblin {
    class Color;
    class Scene;
    class Ray;
    class Renderer {
    public:
        void render(Scene* scene);
        Color Li(Scene* scene, const Ray& ray);
    };
}

#endif //GOBLIN_RENDERER_H