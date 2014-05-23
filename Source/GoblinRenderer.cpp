#include "GoblinRenderer.h"
#include "GoblinRay.h"
#include "GoblinColor.h"
#include "GoblinCamera.h"
#include "GoblinFilm.h"
#include "GoblinUtils.h"

namespace Goblin {

    RenderTask::RenderTask(ImageTile* tile, Renderer* renderer,
        const CameraPtr& camera, const ScenePtr& scene,
        const SampleQuota& sampleQuota, int samplePerPixel): 
        mTile(tile), mRenderer(renderer),
        mCamera(camera), mScene(scene),
        mSampleQuota(sampleQuota), 
        mSamplePerPixel(samplePerPixel) {
        mRNG = new RNG();
    }

    RenderTask::~RenderTask() {
        if(mRNG) {
            delete mRNG;
            mRNG = NULL;
        }
    }

    void RenderTask::run() {
        int xStart, xEnd, yStart, yEnd;
        mTile->getSampleRange(&xStart, &xEnd, &yStart, &yEnd);

        std::cout << "tiles " <<
            " sampleXStart " << xStart << " sampleYStart " << yStart <<
            " sampleXEnd " << xEnd << " sampleYEnd " << yEnd <<
                std::endl;

        Sampler sampler(xStart, xEnd, yStart, yEnd, 
            mSamplePerPixel, mSampleQuota, mRNG);
        int batchAmount = sampler.maxSamplesPerRequest();
        Sample* samples = sampler.allocateSampleBuffer(batchAmount);
        int sampleNum = 0;
        while((sampleNum = sampler.requestSamples(samples)) > 0) {
            for(int s = 0; s < sampleNum; ++s) {
                Ray ray;
                float w = mCamera->generateRay(samples[s], &ray);
                Color L = w * mRenderer->Li(mScene, ray, samples[s], *mRNG);
                mTile->addSample(samples[s], L);
            }
        }
        delete [] samples;
    }

    Renderer::Renderer(const RenderSetting& setting):
        mLightSampleIndexes(NULL), mBSDFSampleIndexes(NULL),
        mPickLightSampleIndexes(NULL),
        mPowerDistribution(NULL), 
        mSetting(setting) {
    }

    Renderer::~Renderer() {
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
    }

    void Renderer::render(const ScenePtr& scene) {
        const CameraPtr camera = scene->getCamera();
        Film* film = camera->getFilm();
        SampleQuota sampleQuota;
        querySampleQuota(scene, &sampleQuota);
        int samplePerPixel = mSetting.samplePerPixel;

        vector<ImageTile*>& tiles = film->getTiles();
        vector<Task*> renderTasks;
        for(size_t i = 0; i < tiles.size(); ++i) {
            renderTasks.push_back(new RenderTask(tiles[i], this, 
                camera, scene, sampleQuota, samplePerPixel));
        }
        
        ThreadPool threadPool;
        threadPool.enqueue(renderTasks);
        threadPool.waitForAll();
        //clean up
        for(size_t i = 0; i < renderTasks.size(); ++i) {
            delete renderTasks[i];
        }
        renderTasks.clear();
        film->writeImage();
    }

    Color Renderer::singleSampleLd(const ScenePtr& scene, const Ray& ray,
        float epsilon, const Intersection& intersection,
        const Sample& sample, 
        const LightSample& lightSample,
        const BSDFSample& bsdfSample,
        float pickLightSample,
        BSDFType type) const {

        const vector<Light*>& lights = scene->getLights();
        if(lights.size() == 0) {
            return Color::Black;
        }
        float pdf;
        int lightIndex = 
            mPowerDistribution->sampleDiscrete(pickLightSample, &pdf);
        const Light* light = lights[lightIndex];
        Color Ld = estimateLd(scene, -ray.d, epsilon, intersection,
            light, lightSample, bsdfSample, type) / pdf;
        return Ld;
    }

    Color Renderer::multiSampleLd(const ScenePtr& scene, const Ray& ray,
        float epsilon, const Intersection& intersection,
        const Sample& sample, const RNG& rng,
        LightSampleIndex* lightSampleIndexes,
        BSDFSampleIndex* bsdfSampleIndexes,
        BSDFType type) const {
        Color totalLd = Color::Black;
        const vector<Light*>& lights = scene->getLights();
        for(size_t i = 0; i < lights.size(); ++i) {
            Color Ld = Color::Black;
            uint32_t samplesNum = lightSampleIndexes[i].samplesNum;
            for(size_t n = 0; n < samplesNum; ++n) {
                const Light* light = lights[i];
                LightSample ls(rng);
                BSDFSample bs(rng);
                if(lightSampleIndexes != NULL && bsdfSampleIndexes != NULL) {
                    ls = LightSample(sample, lightSampleIndexes[i], n);
                    bs = BSDFSample(sample, bsdfSampleIndexes[i], n);
                }
                Ld +=  estimateLd(scene, -ray.d, epsilon, intersection,
                    light, ls, bs, type);
            }
            Ld /= static_cast<float>(samplesNum);
            totalLd += Ld;
        }
        return totalLd;
    }

    Color Renderer::estimateLd(const ScenePtr& scene, const Vector3& wo,
        float epsilon, const Intersection& intersection, const Light* light, 
        const LightSample& ls, const BSDFSample& bs, BSDFType type) const {

        Color Ld = Color::Black;
        const MaterialPtr& material = 
            intersection.primitive->getMaterial();
        const Fragment& fragment = intersection.fragment;
        Vector3 wi;
        Vector3 p = fragment.getPosition();
        Vector3 n = fragment.getNormal();
        float lightPdf, bsdfPdf;
        Ray shadowRay;
        // MIS for lighting part
        Color L = light->sampleL(p, epsilon, ls, &wi, &lightPdf, &shadowRay);
        if(L != Color::Black && lightPdf > 0.0f) {
            Color f = material->bsdf(fragment, wo, wi);
            if(f != Color::Black && !scene->intersect(shadowRay)) {
                // we don't do MIS for delta distribution light
                // since there is only one sample need for it
                if(light->isDelta()) {
                    return f * L * absdot(n, wi) / lightPdf;
                } else {
                    bsdfPdf = material->pdf(fragment, wo, wi);
                    float lWeight = powerHeuristic(1, lightPdf, 1, bsdfPdf);
                    Ld += f * L * absdot(n, wi) * lWeight / lightPdf;
                }
            }
        }

        // MIS for bsdf part
        BSDFType sampledType;
        Color f = material->sampleBSDF(fragment, wo, bs, 
            &wi, &bsdfPdf, type, &sampledType);
        if(f != Color::Black && bsdfPdf > 0.0f) {
            // calculate the misWeight if it's not a specular material
            // otherwise we should got 0 Ld from light sample earlier,
            // and count on this part for all the Ld contribution
            float fWeight = 1.0f;
            if(!(sampledType & BSDFSpecular)) {
                lightPdf = light->pdf(p, wi);
                if(lightPdf == 0.0f) {
                    return Ld;
                }
                fWeight = powerHeuristic(1, bsdfPdf, 1, lightPdf);
            }
            Intersection lightIntersect;
            float lightEpsilon;
            Ray r(fragment.getPosition(), wi, epsilon);
            if(scene->intersect(r, &lightEpsilon, &lightIntersect)) {
                if(lightIntersect.primitive->getAreaLight() == light) {
                    Color Li = lightIntersect.Le(-wi);
                    if(Li != Color::Black) {
                        Ld += f * Li * absdot(wi, n) * fWeight / bsdfPdf;
                    }
                }
            }
        }

        return Ld;
    }

    Color Renderer::specularReflect(const ScenePtr& scene, const Ray& ray, 
        float epsilon, const Intersection& intersection,
        const Sample& sample, const RNG& rng) const {
        Color L(Color::Black);
        const Vector3& n = intersection.fragment.getNormal();
        const Vector3& p = intersection.fragment.getPosition();
        const MaterialPtr& material = 
            intersection.primitive->getMaterial();
        Vector3 wo = -ray.d;
        Vector3 wi;
        float pdf;
        // fill in a random BSDFSample for api request, specular actually
        // don't need to do any monte carlo sampling(only one possible out dir)
        Color f = material->sampleBSDF(intersection.fragment, 
            wo, BSDFSample(rng), &wi, &pdf, 
            BSDFType(BSDFSpecular | BSDFReflection));
        if(f != Color::Black && absdot(wi, n) != 0.0f) {
            Ray reflectiveRay(p, wi, epsilon);
            reflectiveRay.depth = ray.depth + 1;
            Color Lr = Li(scene, reflectiveRay, sample, rng);
            L += f * Lr * absdot(wi, n) / pdf;
        }
        return L;
    }

    Color Renderer::specularRefract(const ScenePtr& scene, const Ray& ray, 
        float epsilon, const Intersection& intersection,
        const Sample& sample, const RNG& rng) const {
        Color L(Color::Black);
        const Vector3& n = intersection.fragment.getNormal();
        const Vector3& p = intersection.fragment.getPosition();
        const MaterialPtr& material = 
            intersection.primitive->getMaterial();
        Vector3 wo = -ray.d;
        Vector3 wi;
        float pdf;
        // fill in a random BSDFSample for api request, specular actually
        // don't need to do any monte carlo sampling(only one possible out dir)
        Color f = material->sampleBSDF(intersection.fragment, 
            wo, BSDFSample(rng), &wi, &pdf, 
            BSDFType(BSDFSpecular | BSDFTransmission));
        if(f != Color::Black && absdot(wi, n) != 0.0f) {
            Ray refractiveRay(p, wi, epsilon);
            refractiveRay.depth = ray.depth + 1;
            Color Lr = Li(scene, refractiveRay, sample, rng);
            L += f * Lr * absdot(wi, n) / pdf;
        }
        return L;
    }
}
