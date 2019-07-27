#include "GoblinWhitted.h"
#include "GoblinRay.h"

namespace Goblin {

WhittedRenderer::WhittedRenderer(int samplePerPixel, int threadNum, 
    int maxRayDepth, int bssrdfSampleNum): 
    Renderer(samplePerPixel, threadNum),
    mMaxRayDepth(maxRayDepth),
    mBssrdfSampleNum(bssrdfSampleNum) {}

WhittedRenderer::~WhittedRenderer() {}

Color WhittedRenderer::Li(const ScenePtr& scene, 
    const RayDifferential& ray, 
    const Sample& sample, const RNG& rng,
    RenderingTLS* tls) const {
    Color Li = Color::Black;
    float epsilon;
    Intersection intersection;
    if (scene->intersect(ray, &epsilon, &intersection)) {
        intersection.computeUVDifferential(ray);
        // if intersect an area light
        Li += intersection.Le(-ray.d);
        // subsurface scattering 
        Li += Lsubsurface(scene, intersection, -ray.d, sample, 
            &mBSSRDFSampleIndex, tls);
        // direct light contribution, for specular part we let
        // specularReflect/specularRefract to deal with it
        Li += multiSampleLd(scene, ray, epsilon, intersection, sample, rng,
            mLightSampleIndexes, mBSDFSampleIndexes,
            BSDFType(BSDFAll & ~BSDFSpecular)); 
        // reflection and refraction
        if (ray.depth < mMaxRayDepth) {
            Li += specularReflect(scene, ray, epsilon, intersection, 
                sample, rng);
            Li += specularRefract(scene, ray, epsilon, intersection,
                sample, rng);
        }
    } else {
        // get image based lighting if ray didn't hit anything
        Li += scene->evalEnvironmentLight(ray);
    }
    return Li;        
}

void WhittedRenderer::querySampleQuota(const ScenePtr& scene, 
        SampleQuota* sampleQuota) {
    if (mLightSampleIndexes) {
        delete [] mLightSampleIndexes;
        mLightSampleIndexes = nullptr;
    } 
    if (mBSDFSampleIndexes) {
        delete [] mBSDFSampleIndexes;
        mBSDFSampleIndexes = nullptr;
    }
    if (mPickLightSampleIndexes) {
        delete [] mPickLightSampleIndexes;
        mPickLightSampleIndexes = nullptr;
    }
    const std::vector<Light*>& lights = scene->getLights();
    mLightSampleIndexes = new LightSampleIndex[lights.size()];
    mBSDFSampleIndexes = new BSDFSampleIndex[lights.size()];
    for (size_t i = 0; i < lights.size(); ++i) {
        uint32_t samplesNum = lights[i]->getSamplesNum();
        mLightSampleIndexes[i] = LightSampleIndex(sampleQuota, samplesNum);
        mBSDFSampleIndexes[i] = BSDFSampleIndex(sampleQuota, samplesNum);
    }
    mPickLightSampleIndexes = new SampleIndex[1];
    mPickLightSampleIndexes[0] = sampleQuota->requestOneDQuota(1);
    mBSSRDFSampleIndex = BSSRDFSampleIndex(sampleQuota, mBssrdfSampleNum);
}

Renderer* createWhitted(const ParamSet& params) {
    int samplePerPixel = params.getInt("sample_per_pixel", 1);
    int threadNum = params.getInt("thread_num", getMaxThreadNum());
    int maxRayDepth = params.getInt("max_ray_depth", 5);
    int bssrdfSampleNum = params.getInt("bssrdf_sample_num", 4);
    return new WhittedRenderer(samplePerPixel, threadNum, 
        maxRayDepth, bssrdfSampleNum);
}
}
