#ifndef GOBLIN_SPPM_H
#define GOBLIN_SPPM_H

#include "GoblinRenderer.h"
#include "GoblinUtils.h"

namespace Goblin {

class SpatialHashGrids;
struct PhotonCache;

struct PixelData {
    PixelData(): Phi(0.0f), Mi(0), throughput(0.0f), pathLength(0),
        Ni(0), Ri(sInvalidRadius), Ld(0.0f), Tau(0.0f), pixelIndex(0) {}

    void reset() {
        Phi = Color(0.0f);
        Mi = 0;
        throughput = Color(0.0f);
        pathLength = 0;
    }
    // data that got reset after each photon pass iteration
    // accumulation of fs * Phi_p (photon contribution)
    Color Phi;
    // number of photon hit this pixel in this photon pass
    int Mi;
    // intersected surface during the ray trace pass
    Fragment fragment;
    // intersected surface bsdf (need this to evaluate bsdf during
    // photon tracing stage)
    const Material* material;
    // the out direction for this visible surface
    // (photon pass need this to evaluate Phi_p)
    Vector3 wo;
    // whether this pixel ray land on surface during ray trace pass
    Color throughput;
    // only contribute GIwhen the phton + ray trace path < max path length
    int pathLength;

    // data that stays through the whole sppm iteration cycle
    // Ni+1 = Ni + alpha * Mi
    float Ni;
    // current photon search radius
    // Ri+1 = Ri * sqrt((Ni + alpha * Mi) / (Ni + Mi))
    float Ri;
    // accumulated direct radiance contribution
    Color Ld;
    // Tau_i+1 = (Tau_i + Phi_i) * (Ri+1 / Ri)^2
    Color Tau;

    size_t pixelIndex;

    static const float sInvalidRadius;
};

class SPPM : public Renderer {
public:
    SPPM(int samplePerPixel, int threadNum, int maxPathLength,
        float mInitialRadius);

    ~SPPM();

    Color Li(const ScenePtr& scene, const RayDifferential& ray, 
        const Sample& sample, const RNG& rng,
        RenderingTLS* tls = nullptr) const;

    void render(const ScenePtr& scene);

    void querySampleQuota(const ScenePtr& scene, SampleQuota* sampleQuota);

    void rayTracePass(const ScenePtr& scene, const Sample& sample,
        int pixelX, int pixelY);

    void photonTracePass(const ScenePtr& scene, const Sample& sample,
        std::vector<PhotonCache>& photonCache);

private:
    int mMaxPathLength;
    std::vector<PixelData> mPixelData;
    SpatialHashGrids* mHashGrids;
    float mInitialRadius;
};

Renderer* createSPPM(const ParamSet &params);

}

#endif // GOBLIN_SPPM_H