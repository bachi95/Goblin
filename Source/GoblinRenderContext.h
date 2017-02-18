#ifndef GOBLIN_RENDER_CONTEXT_H
#define GOBLIN_RENDER_CONTEXT_H
#include "GoblinRenderer.h"
#include "GoblinScene.h"

namespace Goblin {
    class RenderContext {
    public:
        RenderContext(RendererPtr renderer, ScenePtr scene):
          mRenderer(renderer), mScene(scene) {}
        void render();

    public:
        RendererPtr mRenderer;
        ScenePtr mScene;
    };

    inline void RenderContext::render() {
        mRenderer->preprocess(mScene);
        mRenderer->render(mScene);
    }
}

#endif //GOBLIN_RENDER_CONTEXT_H