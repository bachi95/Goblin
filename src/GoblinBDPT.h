#ifndef GOBLIN_BDPT_H
#define GOBLIN_BDPT_H

#include "GoblinFactory.h"
#include "GoblinRenderer.h"
#include "GoblinPathVertex.h"
#include "GoblinUtils.h"

namespace Goblin {

    class BDPT : public Renderer {
    public:
        BDPT(int samplePerPixel, int threadNum,
            int maxPathLength = 5, int debugS = -1, int debugT = -1,
            bool debugNoMIS = false);

        ~BDPT();

        Color Li(const ScenePtr& scene, const RayDifferential& ray, 
            const Sample& sample, const RNG& rng,
            RenderingTLS* tls = NULL) const;

        void render(const ScenePtr& scene);
        
        void evalContribution(const ScenePtr& scene,
            const Sample& sample, const RNG& rng,
            std::vector<PathVertex>& lightPath,
            std::vector<PathVertex>& eyePath,
            std::vector<BDPTMISNode>& misNodes,
            ImageTile* tile) const;

        void querySampleQuota(const ScenePtr& scene,
            SampleQuota* sampleQuota);

    private:
        int constructLightPath(const ScenePtr& scene,
            const Sample& sample, const RNG&rng,
            const std::vector<Light*>& lights,
            std::vector<PathVertex>& lightPath) const;

        int constructEyePath(const ScenePtr& scene,
            const Sample& sample, const RNG& rng,
            const CameraPtr& camera,
            std::vector<PathVertex>& eyePath) const;

        Color evalUnweightedContribution(
            const ScenePtr& scene, const CameraPtr& camera,
            const std::vector<PathVertex>& lightPath, int s,
            const std::vector<PathVertex>& eyePath, int t,
            float& G) const;

        // evaluate the geometry term between two PathVertex
        float evalG(const PathVertex& a, const PathVertex& b) const;

        float evalMIS(const ScenePtr& scene, const CameraPtr& camera,
            const std::vector<PathVertex>& lightPath, int s,
            const std::vector<PathVertex>& eyePath, int t,
            const float Gconnect, std::vector<BDPTMISNode>& misNodes) const;

    private:
        uint64_t mTotalSamplesNum;
        int mMaxPathLength;
        BSDFSampleIndex* mLightPathSampleIndexes;
        BSDFSampleIndex* mEyePathSampleIndexes;
        std::vector<float> mPickLightPdf;
        int mDebugS;
        int mDebugT;
        bool mDebugNoMIS;
    };

    class BDPTCreator : public
        Creator<Renderer, const ParamSet&> {
    public:
        Renderer* create(const ParamSet &params) const;
    };
}

#endif // GOBLIN_BDPT_H