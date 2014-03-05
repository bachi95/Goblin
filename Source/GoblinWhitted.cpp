#include "GoblinWhitted.h"
#include "GoblinRay.h"

namespace Goblin {

    WhittedRenderer::WhittedRenderer(const RenderSetting& setting): 
        Renderer(setting) {}

    WhittedRenderer::~WhittedRenderer() {}

    Color WhittedRenderer::Li(const ScenePtr& scene, const Ray& ray, 
        const Sample& sample) const {
        Color Li = Color::Black;
        float epsilon;
        Intersection intersection;
        if(scene->intersect(ray, &epsilon, &intersection)) {
            // if intersect an area light
            Li += intersection.Le(-ray.d);
            // direct light contribution, for specular part we let
            // specularReflect/specularRefract to deal with it
            Li += multiSampleLd(scene, ray, epsilon, intersection, sample,
                mLightSampleIndexes, mBSDFSampleIndexes,
                BSDFType(BSDFAll & ~BSDFSpecular)); 
            // reflection and refraction
            if(ray.depth < mSetting.maxRayDepth) {
                Li += specularReflect(scene, ray, epsilon, intersection, 
                    sample);
                Li += specularRefract(scene, ray, epsilon, intersection,
                    sample);
            }
        }
        return Li;        
    }

    void WhittedRenderer::querySampleQuota(const ScenePtr& scene, 
            Sampler* sampler) {
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
            mLightSampleIndexes[i] = LightSampleIndex(sampler, samplesNum);
            mBSDFSampleIndexes[i] = BSDFSampleIndex(sampler, 
                samplesNum);
        }
        mPickLightSampleIndexes = new SampleIndex[1];
        mPickLightSampleIndexes[0] = sampler->requestOneDQuota(1);
        vector<float> lightPowers;
        for(size_t i = 0; i < lights.size(); ++i) {
            lightPowers.push_back(lights[i]->power(scene).luminance());
        }
        mPowerDistribution = new CDF1D(lightPowers);
    }
}
