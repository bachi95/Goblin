#include "GoblinAO.h"
#include "GoblinRay.h"

namespace Goblin {

AORenderer::AORenderer(int samplePerPixel, int threadNum,
    int aoSampleNum):
    Renderer(samplePerPixel, threadNum),
    mAOSampleNum(aoSampleNum)
{}

Color AORenderer::Li(const ScenePtr& scene,
    const RayDifferential& ray,
    const Sample& sample, const RNG& rng,
    RenderingTLS* tls) const {
    Color Li = Color::Black;
    float epsilon;
    Intersection intersection;
    if (scene->intersect(ray, &epsilon, &intersection)) {
        const Fragment& fragment = intersection.fragment;
        uint32_t samplesNum = mAOSampleIndex.sampleNum;
        uint32_t occludedNum = 0;
        for (uint32_t n = 0; n < samplesNum; ++n) {
            Vector3 sampleDir = uniformSampleHemisphere(
                sample.u2D[mAOSampleIndex.offset][2 * n],
                sample.u2D[mAOSampleIndex.offset][2 * n + 1]);
            Matrix3 shadeToWorld = fragment.getWorldToShade().transpose();
            Vector3 occludeRayDir = shadeToWorld * sampleDir;
            Ray occludeRay(fragment.getPosition(), occludeRayDir, epsilon);
            if (scene->occluded(occludeRay)) {
                occludedNum++;
            }
        }
        Li = Color((float)(samplesNum - occludedNum) / (float)samplesNum);
    }
    return Li;
}

void AORenderer::querySampleQuota(const ScenePtr& scene,
        SampleQuota* sampleQuota) {
    mAOSampleIndex = sampleQuota->requestTwoDQuota(mAOSampleNum);
}

Renderer* createAO(const ParamSet& params) {
    int samplePerPixel = params.getInt("sample_per_pixel", 1);
    int threadNum = params.getInt("thread_num", getMaxThreadNum());
    int aoSampleNum = params.getInt("ao_sample_num", 25);
    return new AORenderer(samplePerPixel, threadNum, aoSampleNum);
}

}