#include "GoblinPathtracer.h"
#include "GoblinRay.h"

namespace Goblin {

    PathTracer::PathTracer(const RenderSetting& setting): 
        Renderer(setting), mPathSampleIndexes(NULL) {}

    PathTracer::~PathTracer() {
        if(mPathSampleIndexes != NULL) {
            delete [] mPathSampleIndexes;
            mPathSampleIndexes = NULL;
        }
    }

    Color PathTracer::Li(const ScenePtr& scene, const Ray& ray, 
        const Sample& sample, const RNG& rng) const {
        Color Li = Color::Black;
        float epsilon;
        Intersection intersection;
        if(!scene->intersect(ray, &epsilon, &intersection)) {
            // get image based lighting if ray didn't hit anything
            const vector<Light*>& lights = scene->getLights();
            for(size_t i = 0; i < lights.size(); ++i) {
                Li += lights[i]->Le(ray);
            }
            return Li;
        }
        // if intersect an area light
        Li += intersection.Le(-ray.d);

        Ray currentRay = ray;
        Color throughput(1.0f, 1.0f, 1.0f);
        for(int bounces = 0; bounces < mSetting.maxRayDepth; ++bounces) {
            LightSample ls(sample, mLightSampleIndexes[bounces], 0);
            BSDFSample bs(sample, mBSDFSampleIndexes[bounces], 0);
            float pickSample = 
                sample.u1D[mPickLightSampleIndexes[bounces].offset][0];
            Li += throughput * singleSampleLd(scene, 
                currentRay, epsilon, intersection, sample,
                ls, bs, pickSample);

            const MaterialPtr& material = 
                intersection.primitive->getMaterial();
            Vector3 wo = -currentRay.d;
            Vector3 wi;
            float pdf;
            BSDFSample ps(sample, mPathSampleIndexes[bounces], 0); 
            Color f =  material->sampleBSDF(intersection.fragment, wo, ps,
                &wi, &pdf);
            if( f == Color::Black || pdf == 0.0f) {
                break;
            }
            Vector3 p = intersection.fragment.getPosition();
            Vector3 n = intersection.fragment.getNormal();
            throughput *= f * absdot(wi, n) / pdf;
            currentRay = Ray(p, wi, epsilon);
            if(!scene->intersect(currentRay, &epsilon, &intersection)) {
                break;
            }
        }
        return Li;        
    }

    void PathTracer::querySampleQuota(const ScenePtr& scene, 
            SampleQuota* sampleQuota) {
        if(mLightSampleIndexes) {
            delete [] mLightSampleIndexes;
            mLightSampleIndexes = NULL;
        }
        if(mBSDFSampleIndexes) {
            delete [] mBSDFSampleIndexes;
            mBSDFSampleIndexes = NULL;
        }
        if(mPowerDistribution) {
            delete mPowerDistribution;
            mPowerDistribution = NULL;
        }
        if(mPickLightSampleIndexes) {
            delete [] mPickLightSampleIndexes;
            mPickLightSampleIndexes = NULL;
        }
        if(mPathSampleIndexes) {
            delete [] mPathSampleIndexes;
            mPathSampleIndexes = NULL;
        }

        int bounces = mSetting.maxRayDepth;
        mLightSampleIndexes = new LightSampleIndex[bounces];
        mBSDFSampleIndexes = new BSDFSampleIndex[bounces];
        mPathSampleIndexes = new BSDFSampleIndex[bounces];
        mPickLightSampleIndexes = new SampleIndex[bounces];
        for(int i = 0; i < bounces; ++i) {
            mLightSampleIndexes[i] = LightSampleIndex(sampleQuota, 1);
            mBSDFSampleIndexes[i] = BSDFSampleIndex(sampleQuota, 1);
            mPathSampleIndexes[i] = BSDFSampleIndex(sampleQuota, 1);
            mPickLightSampleIndexes[i] = sampleQuota->requestOneDQuota(1);
        }

        const vector<Light*>& lights = scene->getLights();
        vector<float> lightPowers;
        for(size_t i = 0; i < lights.size(); ++i) {
            lightPowers.push_back(lights[i]->power(scene).luminance());
        }
        mPowerDistribution = new CDF1D(lightPowers);
    }
}
