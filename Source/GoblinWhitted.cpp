#include "GoblinWhitted.h"
#include "GoblinRay.h"

namespace Goblin {

    WhittedRenderer::WhittedRenderer(const RenderSetting& setting): 
        Renderer(setting) {}

    WhittedRenderer::~WhittedRenderer() {}

    Color WhittedRenderer::Li(const ScenePtr& scene, const Ray& ray, 
        const Sample& sample, const RNG& rng) const {
        Color Li = Color::Black;
        float epsilon;
        Intersection intersection;
        if(scene->intersect(ray, &epsilon, &intersection)) {
            // if intersect an area light
            Li += intersection.Le(-ray.d);
            // direct light contribution, for specular part we let
            // specularReflect/specularRefract to deal with it
            Li += multiSampleLd(scene, ray, epsilon, intersection, sample, rng,
                mLightSampleIndexes, mBSDFSampleIndexes,
                BSDFType(BSDFAll & ~BSDFSpecular)); 
            // reflection and refraction
            if(ray.depth < mSetting.maxRayDepth) {
                Li += specularReflect(scene, ray, epsilon, intersection, 
                    sample, rng);
                Li += specularRefract(scene, ray, epsilon, intersection,
                    sample, rng);
            }
        } else {
            // get image based lighting if ray didn't hit anything
            const vector<Light*>& lights = scene->getLights();
            for(size_t i = 0; i < lights.size(); ++i) {
                Li += lights[i]->Le(ray);
            }
            return Li;
        }
        return Li;        
    }

    void WhittedRenderer::querySampleQuota(const ScenePtr& scene, 
            SampleQuota* sampleQuota) {
        if(mLightSampleIndexes) {
            delete [] mLightSampleIndexes;
            mLightSampleIndexes = NULL;
        } 
        if(mBSDFSampleIndexes) {
            delete [] mBSDFSampleIndexes;
            mBSDFSampleIndexes = NULL;
        }
        if(mPickLightSampleIndexes) {
            delete [] mPickLightSampleIndexes;
            mPickLightSampleIndexes = NULL;
        }
        if(mPowerDistribution) {
            delete mPowerDistribution;
            mPowerDistribution = NULL;
        }
        const vector<Light*>& lights = scene->getLights();
        mLightSampleIndexes = new LightSampleIndex[lights.size()];
        mBSDFSampleIndexes = new BSDFSampleIndex[lights.size()];
        for(size_t i = 0; i < lights.size(); ++i) {
            uint32_t samplesNum = lights[i]->getSamplesNum();
            mLightSampleIndexes[i] = LightSampleIndex(sampleQuota, samplesNum);
            mBSDFSampleIndexes[i] = BSDFSampleIndex(sampleQuota, samplesNum);
        }
        mPickLightSampleIndexes = new SampleIndex[1];
        mPickLightSampleIndexes[0] = sampleQuota->requestOneDQuota(1);
        vector<float> lightPowers;
        for(size_t i = 0; i < lights.size(); ++i) {
            lightPowers.push_back(lights[i]->power(scene).luminance());
        }
        mPowerDistribution = new CDF1D(lightPowers);
    }
}
