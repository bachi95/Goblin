#include "GoblinRenderer.h"
#include "GoblinRay.h"
#include "GoblinColor.h"
#include "GoblinCamera.h"
#include "GoblinFilm.h"
#include "GoblinUtils.h"

namespace Goblin {

    Renderer::Renderer(const RenderSetting& setting):
        mLightSampleIndexes(NULL),
        mSamples(NULL), mSampler(NULL), mSetting(setting) {
    }

    Renderer::~Renderer() {
        if(mLightSampleIndexes) {
            delete [] mLightSampleIndexes;
            mLightSampleIndexes = NULL;
        }
        if(mSamples) {
            delete [] mSamples;
            mSamples = NULL;
        }
        if(mSampler) {
            delete mSampler;
            mSampler = NULL;
        }
    }

    void Renderer::querySampleQuota(const ScenePtr& scene, Sampler* sampler) {
        if(mLightSampleIndexes) {
            delete [] mLightSampleIndexes;
            mLightSampleIndexes = NULL;
        } 

        const vector<Light*>& lights = scene->getLights();
        mLightSampleIndexes = new LightSampleIndex[lights.size()];
        for(size_t i = 0; i < lights.size(); ++i) {
            uint32_t samplesNum = lights[i]->getSamplesNum();
            mLightSampleIndexes[i] = LightSampleIndex(sampler, samplesNum);
        }
    }

    void Renderer::render(const ScenePtr& scene) {
        const CameraPtr camera = scene->getCamera();
        Film* film = camera->getFilm();
        int xStart, xEnd, yStart, yEnd;
        film->getSampleRange(&xStart, &xEnd, &yStart, &yEnd);
        if(mSampler != NULL) {
            delete mSampler;
        }
        if(mSamples != NULL) {
            delete [] mSamples;
        }
        int samplePerPixel = mSetting.samplePerPixel;
        mSampler = new Sampler(xStart, xEnd, yStart, yEnd, 
            samplePerPixel);
        // bluh....api dependency on the above sampler
        // TODO wragle out the on APIdependency  
        querySampleQuota(scene, mSampler);
        int batchAmount = mSampler->maxSamplesPerRequest();
        mSamples = mSampler->allocateSampleBuffer(batchAmount);
        // temp progress reporter so that waiting can be not that boring...
        // TODO make this something more elegant
        unsigned long int accumulatedBuffer = 0;
        unsigned long int accumulatedSamples = 0;
        unsigned long int maxTotalSamples = mSampler->maxTotalSamples();
        unsigned long int reportStep = maxTotalSamples / 100;
        string backspace = "\b\b\b\b\b\b\b\b\b\b\b\b\b";

        int sampleNum = 0;
        while((sampleNum = mSampler->requestSamples(mSamples)) > 0) {
            for(int i = 0; i < sampleNum; ++i) {
                Ray ray;
                float w = camera->generateRay(mSamples[i], &ray);
                Color L = w * Li(scene, ray, mSamples[i]);
                film->addSample(mSamples[i], L);
            }

            // print out progress
            accumulatedBuffer += sampleNum;
            if(accumulatedBuffer > reportStep) {
                accumulatedSamples += accumulatedBuffer;
                accumulatedBuffer = 0;
                std::cout << backspace;
                std::cout << "progress %";
                std::cout << 100 * accumulatedSamples / maxTotalSamples;
                std::cout.flush();
            }
        }
        std::cout << backspace;

        //for now leave this here as debug purpose......
        //Ray ray(Vector3(0, 0, 0.0), Vector3(0, 0, 1), 0); 
        //Color L = Li(scene, ray);

        film->writeImage();
    }

    Color Renderer::Li(const ScenePtr& scene, const Ray& ray, 
        const Sample& sample) const {
        Color Li = Color::Black;
        float epsilon;
        Intersection intersection;
        if(scene->intersect(ray, &epsilon, &intersection)) {
            // if intersect an area light
            Li += intersection.Le(-ray.d);
            // direct light contribution
            Li += directLighting(scene, ray, epsilon, intersection, sample); 
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

    Color Renderer::directLighting(const ScenePtr& scene, const Ray& ray,
        float epsilon, const Intersection& intersection,
        const Sample& sample) const {

        Color totalLd = Color::Black;
        const vector<Light*>& lights = scene->getLights();
        for(size_t i = 0; i < lights.size(); ++i) {
            Color Ld = Color::Black;
            uint32_t samplesNum = mLightSampleIndexes[i].samplesNum;
            for(size_t n = 0; n < samplesNum; ++n) {
                const Light* light = lights[i];
                LightSample ls(sample, mLightSampleIndexes[i], n);
                Ld +=  estimateLd(scene, ray, epsilon, intersection,
                    light, ls);
            }
            Ld /= static_cast<float>(samplesNum);
            totalLd += Ld;
        }
        return totalLd;
    }


    Color Renderer::estimateLd(const ScenePtr& scene, const Ray& ray,
        float epsilon, const Intersection& intersection, 
        const Light* light, const LightSample& ls) const {

        Color Ld = Color::Black;
        const MaterialPtr& material = 
            intersection.primitive->getMaterial();
        const Fragment& fragment = intersection.fragment;
        Vector3 wi;
        float pdf;
        Ray shadowRay;
        Color L = light->sampleL(fragment.position, epsilon, 
            ls, &wi, &pdf, &shadowRay);

        if(L != Color::Black && pdf > 0.0f) {
            Vector3 wo = -ray.d;
            Color f = material->bsdf(fragment, wo, wi);
            if(f != Color::Black && !scene->intersect(shadowRay)) {
                Ld += f * L * absdot(fragment.normal, wi) / pdf;
            }
        }

        return Ld;
    }

    Color Renderer::specularReflect(const ScenePtr& scene, const Ray& ray, 
        float epsilon, const Intersection& intersection,
        const Sample& sample) const {
        Color L(Color::Black);
        const Vector3& n = intersection.fragment.normal;
        const Vector3& p = intersection.fragment.position;
        const MaterialPtr& material = 
            intersection.primitive->getMaterial();
        Vector3 wo = -ray.d;
        Vector3 wi;
        Color f = material->sampleBSDF(intersection.fragment, 
            wo, &wi, BSDFType(BSDFSpecular | BSDFReflection));
        if(f != Color::Black && absdot(wi, n) != 0.0f) {
            Ray reflectiveRay(p, wi, epsilon);
            reflectiveRay.depth = ray.depth + 1;
            Color Lr = Li(scene, reflectiveRay, sample);
            L += f * Lr * absdot(wi, n);
        }
        return L;
    }

    Color Renderer::specularRefract(const ScenePtr& scene, const Ray& ray, 
        float epsilon, const Intersection& intersection,
        const Sample& sample) const {
        Color L(Color::Black);
        const Vector3& n = intersection.fragment.normal;
        const Vector3& p = intersection.fragment.position;
        const MaterialPtr& material = 
            intersection.primitive->getMaterial();
        Vector3 wo = -ray.d;
        Vector3 wi;
        Color f = material->sampleBSDF(intersection.fragment, 
            wo, &wi, BSDFType(BSDFSpecular | BSDFTransmission));
        if(f != Color::Black && absdot(wi, n) != 0.0f) {
            Ray refractiveRay(p, wi, epsilon);
            refractiveRay.depth = ray.depth + 1;
            Color Lr = Li(scene, refractiveRay, sample);
            L += f * Lr * absdot(wi, n);
        }
        return L;
    }

}
