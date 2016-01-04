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
            throughput(0.0f),
            light(NULL), camera(NULL), material(NULL) {}

        PathVertex(const Color& t, const Vector3& p, const Vector3& n,
            const Light* l): throughput(t),
            fragment(p, n, Vector2::Zero, Vector3::Zero, Vector3::Zero),
            light(l), camera(NULL), material(NULL) {}

        PathVertex(const Color& t, const Vector3& p, const Vector3& n,
            const Camera* c): throughput(t),
            fragment(p, n, Vector2::Zero, Vector3::Zero, Vector3::Zero),
            light(NULL), camera(c), material(NULL) {}

        PathVertex(const Color& t, const Fragment& f,
            const Material* m): throughput(t), fragment(f),
            light(NULL), camera(NULL), material(m) {}

        Color throughput;
        Fragment fragment;
        const Light* light;
        const Camera* camera;
        const Material* material;
    };

    class LightTracer : public Renderer {
    public:
        LightTracer(int samplePerPixel, int threadNum,
            int maxPathLength = 5);

        ~LightTracer();

        Color Li(const ScenePtr& scene, const RayDifferential& ray, 
            const Sample& sample, const RNG& rng,
            WorldDebugData* debugData = NULL) const;

        void render(const ScenePtr& scene);
        
        // t = 1 strategy
        // random walk a particle path from light source and connect to
        // one camera sample, kinda like the adjoint version of common
        // path tracing technique
        void splatFilmT1(const ScenePtr& scene, const Sample& sample,
            const RNG& rng, std::vector<PathVertex>& pathVertices,
            ImageTile* tile) const;
        // t = 0 strategy
        // random walk a particle path from light source and only contribute
        // to film when the last intersection hit the camera lens. This is
        // probably the worst strategy to approximate engergy integration.
        // Still implement it to verify we can converge to the same result
        // compare to other strategy
        void splatFilmT0(const ScenePtr& scene, const Sample& sample,
            const RNG& rng, std::vector<PathVertex>& pathVertices,
            ImageTile* tile) const;

        // s = 1 strategy
        // random walk a particle path from camera and and connect to one
        // light source sample, this is under the hood how path tracing is
        // implemented.
        void splatFilmS1(const ScenePtr& scene, const Sample& sample,
            const RNG& rng, std::vector<PathVertex>& pathVertices,
            ImageTile* tile) const;

        void querySampleQuota(const ScenePtr& scene,
            SampleQuota* sampleQuota);
    private:
        uint64_t mTotalSamplesNum;
        int mMaxPathLength;
    };

    class LightTracerCreator : public
        Creator<Renderer, const ParamSet&> {
    public:
        Renderer* create(const ParamSet &params) const;
    };
}

#endif // GOBLIN_LIGHTTRACER_H