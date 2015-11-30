#ifndef GOBLIN_LIGHTTRACER_H
#define GOBLIN_LIGHTTRACER_H

#include "GoblinFactory.h"
#include "GoblinRenderer.h"
#include "GoblinLight.h"
#include "GoblinMaterial.h"
#include "GoblinCamera.h"
#include "GoblinUtils.h"

namespace Goblin {
    class PathVertex{
    public:
        PathVertex():
            throughput(0.0f), position(Vector3::Zero), normal(Vector3::Zero),
            light(NULL), camera(NULL), material(NULL) {}

        PathVertex(const Color& t, const Vector3& p, const Vector3& n,
            const Light* l): throughput(t), position(p), normal(n),
            light(l), camera(NULL), material(NULL) {}

        PathVertex(const Color& t, const Vector3& p, const Vector3& n,
            const Camera* c): throughput(t), position(p), normal(n),
            light(NULL), camera(c), material(NULL) {}

        PathVertex(const Color& t, const Vector3& p, const Vector3& n,
            const Material* m): throughput(t), position(p), normal(n),
            light(NULL), camera(NULL), material(m) {}

        Color throughput;
        Vector3 position;
        Vector3 normal;
        const Light* light;
        const Camera* camera;
        const Material* material;
    };

    class LightTracer : public Renderer {
    public:
        LightTracer(int samplePerPixel = 1,
            int maxPathLength = 5);

        ~LightTracer();

        Color Li(const ScenePtr& scene, const RayDifferential& ray, 
            const Sample& sample, const RNG& rng,
            WorldDebugData* debugData = NULL) const;

        void render(const ScenePtr& scene);
    private:
        
        void splatFilm(const ScenePtr& scene, const Sample& sample,
            const RNG& rng);

        void querySampleQuota(const ScenePtr& scene,
            SampleQuota* sampleQuota);
    private:
        uint64_t mTotalSamplesNum;
        int mMaxPathLength;
        std::vector<PathVertex> mPathVertices;
    };

    class LightTracerCreator : public
        Creator<Renderer, const ParamSet&> {
    public:
        Renderer* create(const ParamSet &params) const;
    };
}

#endif // GOBLIN_LIGHTTRACER_H