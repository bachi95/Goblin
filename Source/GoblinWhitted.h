#ifndef GOBLIN_WHITTED_H
#define GOBLIN_WHITTED_H

#include "GoblinRenderer.h"

namespace Goblin {
    class WhittedRenderer : public Renderer {
    public:
        WhittedRenderer(const RenderSetting& setting);
        ~WhittedRenderer();
        Color Li(const ScenePtr& scene, const Ray& ray, 
            const Sample& sample) const;
    private:
        void querySampleQuota(const ScenePtr& scene, 
            Sampler* sampler);
    };
}

#endif //GOBLIN_WHITTED_H
