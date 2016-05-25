#include "GoblinRenderer.h"
#include "GoblinRay.h"
#include "GoblinColor.h"
#include "GoblinCamera.h"
#include "GoblinFilm.h"
#include "GoblinUtils.h"
#include "GoblinVolume.h"

namespace Goblin {

    RenderTask::RenderTask(Renderer* renderer, const CameraPtr& camera,
        const ScenePtr& scene, const SampleRange& sampleRange,
        const SampleQuota& sampleQuota, int samplePerPixel,
        RenderProgress* renderProgress): 
        mRenderer(renderer), mCamera(camera), mScene(scene),
        mSampleRange(sampleRange), mSampleQuota(sampleQuota), 
        mSamplePerPixel(samplePerPixel),
        mRenderProgress(renderProgress) {
        mRNG = new RNG();
    }

    RenderTask::~RenderTask() {
        if(mRNG) {
            delete mRNG;
            mRNG = NULL;
        }
    }

    void RenderTask::run(TLSPtr& tls) {
        RenderingTLS* renderingTLS =
            static_cast<RenderingTLS*>(tls.get());
        ImageTile* tile = renderingTLS->getTile();

        Sampler sampler(mSampleRange, mSamplePerPixel, mSampleQuota, mRNG);
        int batchAmount = sampler.maxSamplesPerRequest();
        Sample* samples = sampler.allocateSampleBuffer(batchAmount);
        int sampleNum = 0;
        while((sampleNum = sampler.requestSamples(samples)) > 0) {
            for(int s = 0; s < sampleNum; ++s) {
                RayDifferential ray;
                float w = mCamera->generateRay(samples[s], &ray);
                Color L = mRenderer->Li(mScene, ray, samples[s], 
                    *mRNG, renderingTLS);
                Color tr = mRenderer->transmittance(mScene, ray);
                Color Lv = mRenderer->Lv(mScene, ray, *mRNG);
                tile->addSample(samples[s].imageX, samples[s].imageY,
                    w * (tr * L + Lv));
            }
        }
        delete [] samples;
        mRenderProgress->update();
    }

    RenderProgress::RenderProgress(int taskNum): 
        mFinishedNum(0), mTasksNum(taskNum) {
    }

    void RenderProgress::reset() {
        boost::lock_guard<boost::mutex> lk(mUpdateMutex);
        mFinishedNum = 0;
    }

    void RenderProgress::update() {
        boost::lock_guard<boost::mutex> lk(mUpdateMutex);
        mFinishedNum++;
        std::cout.precision(3);
        std::cout << "\rProgress: %" << 
            (float)mFinishedNum / (float)mTasksNum * 100.0f <<
            "                     ";

        std::cout.flush();
        if(mFinishedNum == mTasksNum) {
            std::cout << "\rRender Complete!         " << std::endl;
            std::cout.flush();
        }
    }

    Renderer::Renderer(int samplePerPixel, int threadNum):
        mLightSampleIndexes(NULL), mBSDFSampleIndexes(NULL),
        mPickLightSampleIndexes(NULL),
        mPowerDistribution(NULL), 
        mSamplePerPixel(samplePerPixel),
        mThreadNum(threadNum) {}

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

        vector<SampleRange> sampleRanges;
        getSampleRanges(film, sampleRanges);
        vector<Task*> renderTasks;
        RenderProgress progress(sampleRanges.size());
        for(size_t i = 0; i < sampleRanges.size(); ++i) {
            renderTasks.push_back(new RenderTask(this, 
                camera, scene, sampleRanges[i], sampleQuota, mSamplePerPixel,
                &progress));
        }
        
        RenderingTLSManager tlsManager(film);
        ThreadPool threadPool(mThreadNum, &tlsManager);
        threadPool.enqueue(renderTasks);
        threadPool.waitForAll();
        //clean up
        for(size_t i = 0; i < renderTasks.size(); ++i) {
            delete renderTasks[i];
        }
        renderTasks.clear();
        drawDebugData(tlsManager.getDebugData(), camera);
        film->writeImage();
    }

    Color Renderer::LbssrdfSingle(const ScenePtr& scene,
        const Fragment& fragment, const BSSRDF* bssrdf, const Vector3& wo,
        const Sample& sample, 
        const BSSRDFSampleIndex* bssrdfSampleIndex,
        RenderingTLS* tls) const {
        const Vector3& pwo = fragment.getPosition();
        const Vector3& no = fragment.getNormal();
        float coso = absdot(wo, fragment.getNormal());
        float eta = bssrdf->getEta();
        float Ft = 1.0f - Material::fresnelDieletric(coso, 1.0f, eta);
        Color scatter = bssrdf->getScatter(fragment);
        Color sigmaT = bssrdf->getAttenuation(fragment);
        float falloff = sigmaT.luminance();
        Vector3 woRefract = Goblin::specularRefract(wo, no, 1.0f, eta);
        const vector<Light*>& lights = scene->getLights();
        Color Lsinglescatter(0.0f);
        // single scattering part
        for(uint32_t i = 0; i < bssrdfSampleIndex->samplesNum; ++i) {
            const BSSRDFSample bssrdfSample(sample, *bssrdfSampleIndex, i);
            // sample a distance with exponential falloff
            float d = exponentialSample(bssrdfSample.uSingleScatter, falloff);
            Vector3 pSample = pwo + d * woRefract;
            float samplePdf = exponentialPdf(d, falloff);
            // sample a light from pSample
            float pickLightPdf;
            int lightIndex = mPowerDistribution->sampleDiscrete(
                bssrdfSample.uPickLight, &pickLightPdf);
            const Light* light = lights[lightIndex];
            Vector3 wi;
            float epsilon, lightPdf;
            Ray shadowRay;
            Color L = light->sampleL(pSample, 1e-5f, bssrdfSample.ls, 
                &wi, &lightPdf, &shadowRay);
            if(L == Color::Black || lightPdf == 0.0f) {
                continue;
            }
            // preserve the maxt since the next intersection test will modify
            // it to the closest intersection, but we still need this maxt to 
            // cast shadow ray from that intersection to light
            float maxt = shadowRay.maxt;
            Intersection wiIntersect;
            // TODO we can use a material indexed lookup table to query
            // surface with this particular BSSRDF instead of doing a whole
            // scene intersaction test 
            if(scene->intersect(shadowRay, &epsilon, &wiIntersect)) {
                // if not, the pSample is out of the BSSRDF geometry already
                if(wiIntersect.getMaterial()->getBSSRDF() == bssrdf) {
                    // update shadow ray to start from pwi
                    const Fragment& fwi = wiIntersect.fragment;
                    const Vector3& pwi = fwi.getPosition();
                    const Vector3& ni = fwi.getNormal();
                    shadowRay.mint = shadowRay.maxt + epsilon;
                    shadowRay.maxt = maxt;
                    if(!scene->intersect(shadowRay)) {
                        float p = bssrdf->phase(wi, woRefract); 
                        float cosi = absdot(ni, wi);
                        float Fti = 1.0f - 
                            Material::fresnelDieletric(cosi, 1.0f, eta);
                        Color sigmaTi = bssrdf->getAttenuation(fwi);
                        float G = absdot(ni, woRefract) / cosi;
                        Color sigmaTC =  sigmaT + G * sigmaTi;
                        float di = length(pwi - pSample);
                        float et = 1.0f / eta;
                        float diPrime = di * absdot(wi, ni) / 
                            sqrt(1.0f - et * et * (1.0f - cosi * cosi));
                        Lsinglescatter += (Ft * Fti * p * scatter / sigmaTC) *
                            expColor(-diPrime * sigmaTi) * 
                            expColor(-d * sigmaT) * L /
                             (lightPdf * pickLightPdf * samplePdf); 
                    }
                }
            } 
        } 
        Lsinglescatter /= (float)bssrdfSampleIndex->samplesNum; 
        return Lsinglescatter;
    }

    Color Renderer::LbssrdfDiffusion(const ScenePtr& scene,
        const Fragment& fragment, const BSSRDF* bssrdf, const Vector3& wo,
        const Sample& sample, 
        const BSSRDFSampleIndex* bssrdfSampleIndex,
        RenderingTLS* tls) const {
        const Vector3& pwo = fragment.getPosition();
        float coso = absdot(wo, fragment.getNormal());
        float eta = bssrdf->getEta();
        float Ft = 1.0f - Material::fresnelDieletric(coso, 1.0f, eta);
        float sigmaTr = bssrdf->getSigmaTr(fragment).luminance();
        const vector<Light*>& lights = scene->getLights();
        // figure out the sample probe radius, we ignore the integration
        // of area with pdf too small compare to center, yes...this introduces
        // bias...but can limit the probe ray as short as possible
        // skipRatio = pdf(r) / pdf(0, 0) = exp(-sigmaTr * r^2)
        float skipRatio = 0.01f;
        float Rmax = sqrt(log(skipRatio) / -sigmaTr);
        Color Lmultiscatter(0.0f);
        for(uint32_t i = 0; i < bssrdfSampleIndex->samplesNum; ++i) {
            const BSSRDFSample bssrdfSample(sample, *bssrdfSampleIndex, i);
            // sample a probe ray with gaussian falloff pdf from intersection
            Ray probeRay;
            float discPdf;
            BSSRDFSampleAxis axis = bssrdf->sampleProbeRay(fragment, 
                bssrdfSample, sigmaTr, Rmax, &probeRay, &discPdf);
            Intersection probeIntersect;
            float epsilon;
            // TODO we can use a material indexed lookup table to query
            // surface with this particular BSSRDF instead of doing a whole
            // scene intersaction test 
            if(scene->intersect(probeRay, &epsilon, &probeIntersect)) {
                if(probeIntersect.getMaterial()->getBSSRDF() == bssrdf) {
                    const Fragment& probeFragment = probeIntersect.fragment;
                    const Vector3& pProbe = probeFragment.getPosition();
                    Color Rd = bssrdf->Rd(probeFragment, 
                        squaredLength(pProbe - pwo));
                    // calculate the irradiance on the sample point
                    float pickLightPdf;
                    int lightIndex = mPowerDistribution->sampleDiscrete(
                        bssrdfSample.uPickLight, &pickLightPdf);
                    const Light* light = lights[lightIndex];
                    Vector3 wi;
                    float lightPdf;
                    Ray shadowRay;
                    const Vector3& ni = probeFragment.getNormal();
                    Color L = light->sampleL(pProbe, epsilon, bssrdfSample.ls, 
                        &wi, &lightPdf, &shadowRay);
                    if(L == Color::Black || lightPdf == 0.0f ||
                        scene->intersect(shadowRay)) {
                        continue;
                    }
                    float cosi = absdot(ni, wi);
                    Color irradiance = L * cosi / (lightPdf * pickLightPdf);
                    float Fti = 1.0f - 
                        Material::fresnelDieletric(cosi, 1.0f, eta);
                    // evaluate the MIS weight
                    float pdf = discPdf * absdot(probeRay.d, ni);
                    float w = bssrdf->MISWeight(fragment, probeFragment, axis, 
                        pdf, sigmaTr, Rmax);
                    Lmultiscatter += 
                        (w * INV_PI * Ft  * Fti * Rd * irradiance) / pdf;
                }
            }

        }
        Lmultiscatter /= (float)bssrdfSampleIndex->samplesNum;
        return Lmultiscatter;
    }

    Color Renderer::Lsubsurface(const ScenePtr& scene,
        const Intersection& intersection, const Vector3& wo,
        const Sample& sample, const BSSRDFSampleIndex* bssrdfSampleIndex,
        RenderingTLS* tls) const {
        const MaterialPtr& material = 
            intersection.primitive->getMaterial();
        const BSSRDF* bssrdf = material->getBSSRDF();
        const vector<Light*>& lights = scene->getLights();
        if(bssrdf == NULL || lights.size() == 0) {
            return Color(0.0f);
        }
        const Fragment& fragment = intersection.fragment;
        Color Lsinglescatter = LbssrdfSingle(scene, fragment, bssrdf, 
            wo, sample, bssrdfSampleIndex, tls);
        // multiple scattering part with diffusion approximation
        Color Lmultiscatter = LbssrdfDiffusion(scene, fragment, bssrdf, 
            wo, sample, bssrdfSampleIndex, tls);
        return Lsinglescatter + Lmultiscatter;
    }

    Color Renderer::Lv(const ScenePtr& scene, const Ray& ray, 
        const RNG& rng) const {
        const VolumeRegion* volume = scene->getVolumeRegion();
        float tMin, tMax;
        if(!volume || !volume->intersect(ray, &tMin, &tMax)) {
            return Color::Black;
        }
        float stepSize = volume->getSampleStepSize();
        Vector3 pPrevious = ray(tMin);
        float tCurrent = tMin + stepSize * rng.randomFloat();
        Vector3 pCurrent = ray(tCurrent);

        Color Lv(0.0f);
        Color transmittance(1.0f);
        while(tCurrent <= tMax) {
            Ray rSegment(pPrevious, pCurrent - pPrevious, 0.0f, 1.0f);
            Color trSegment = volume->transmittance(rSegment);
            transmittance *= trSegment;
            // the emission part
            Lv += transmittance * volume->getEmission(pCurrent);
            // sample light for in scattring part
            Color scatter = volume->getScatter(pCurrent);
            float pickLightSample = rng.randomFloat();
            float pickLightPdf;
            int lightIndex = mPowerDistribution->sampleDiscrete(
                pickLightSample, &pickLightPdf);
            const vector<Light*>& lights = scene->getLights();
            if(lights.size() > 0) {
                const Light* light = lights[lightIndex];
                Ray shadowRay;
                Vector3 wi;
                float lightPdf;
                LightSample ls(rng);
                Color L = light->sampleL(pCurrent, 0.0f, ls, 
                    &wi, &lightPdf, &shadowRay);
                if(L != Color::Black && lightPdf > 0.0f) {
                    if(!scene->intersect(shadowRay)) {
                        Color Ld = volume->transmittance(shadowRay) * L / 
                            (pickLightPdf * lightPdf);
                        float phase = volume->phase(pCurrent, ray.d, wi);
                        Lv += transmittance * scatter * phase * Ld;
                    }
                }
            }
            // advance to the next sample segment
            tCurrent += stepSize;
            pPrevious = pCurrent;
            pCurrent = ray(tCurrent);
        }
        /*
         * the monte carlo estimator for integrate source term
         * from tMin to tMax is (1 / N) * sum(source_term(pi), 1, N) / pdf(pi)
         * where pdf(pi) is 1 / (tMax - tMin) and stepSize = (tMax - tMin) / N
         * which give us the following sweet result
         */
        return stepSize * Lv;

    }

    Color Renderer::transmittance(const ScenePtr& scene, 
        const Ray& ray) const {
        const VolumeRegion* volume = scene->getVolumeRegion();
        if(volume == NULL) {
            return Color(1.0f);
        }
        return volume->transmittance(ray);
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
                Ld += estimateLd(scene, -ray.d, epsilon, intersection,
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
            } else {
                // the radiance contribution from IBL
                Ld += f * light->Le(r, bsdfPdf, sampledType) * 
                    fWeight / bsdfPdf;
            }
        }
        return Ld;
    }

    Color Renderer::singleSampleIrradiance(const ScenePtr& scene,
        float epsilon, const Intersection& intersection,
        const LightSample& lightSample, float pickLightSample) const {
        const vector<Light*>& lights = scene->getLights();
        if(lights.size() == 0) {
            return Color::Black;
        }
        float pdf;
        int lightIndex = 
            mPowerDistribution->sampleDiscrete(pickLightSample, &pdf);
        const Light* light = lights[lightIndex];
        Color irradiance = estimateIrradiance(scene, epsilon, intersection,
            light, lightSample) / pdf;
        return irradiance;
    }

    Color Renderer::estimateIrradiance(const ScenePtr& scene,
        float epsilon, const Intersection& intersection, 
        const Light* light, const LightSample& ls) const {
        Color irradiance(0.0f);
        const Fragment& fragment = intersection.fragment;
        Vector3 wi;
        const Vector3& p = fragment.getPosition();
        const Vector3& n = fragment.getNormal();
        float pdf;
        Ray shadowRay;
        Color L = light->sampleL(p, epsilon, ls, &wi, &pdf, &shadowRay);
        if(L != Color::Black && pdf > 0.0f) {
            if(!scene->intersect(shadowRay)) {
                irradiance += L * absdot(n, wi) / pdf;
            }
        }
        return irradiance;
    }

    Color Renderer::specularReflect(const ScenePtr& scene, 
        const RayDifferential& ray, 
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
            RayDifferential reflectiveRay(p, wi, epsilon);
            reflectiveRay.depth = ray.depth + 1;
            Color Lr = Li(scene, reflectiveRay, sample, rng);
            L += f * Lr * absdot(wi, n) / pdf;
        }
        return L;
    }

    Color Renderer::specularRefract(const ScenePtr& scene, 
        const RayDifferential& ray, 
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
            RayDifferential refractiveRay(p, wi, epsilon);
            refractiveRay.depth = ray.depth + 1;
            Color Lr = Li(scene, refractiveRay, sample, rng);
            L += f * Lr * absdot(wi, n) / pdf;
        }
        return L;
    }

    void Renderer::getSampleRanges(const Film* film,
        vector<SampleRange>& sampleRanges) const {
        SampleRange fullRange;
        film->getSampleRange(fullRange);
        int step = 8;
        for (int y = fullRange.yStart; y < fullRange.yEnd; y += step) {
            for (int x = fullRange.xStart; x < fullRange.xEnd; x += step) {
                SampleRange subSampleRange;
                subSampleRange.xStart = x;
                subSampleRange.xEnd = min(x + step, fullRange.xEnd);
                subSampleRange.yStart = y;
                subSampleRange.yEnd = min(y + step, fullRange.yEnd);
                sampleRanges.push_back(subSampleRange);
            }
        }
    }

    void Renderer::drawDebugData(const DebugData& debugData,
        const CameraPtr& camera) const {
        Film* film = camera->getFilm();
        // if there is any debug data, transform it to screen space
        // and inject to image tile
        const vector<pair<Ray, Color> >& debugRays = debugData.getRays();
        for(size_t i = 0; i < debugRays.size(); ++i) {
            const Ray& r = debugRays[i].first;
            Vector3 sWorld = r.o;
            Vector3 eWorld = r(r.maxt);
            Vector3 sScreen = camera->worldToScreen(sWorld);
            Vector3 eScreen = camera->worldToScreen(eWorld);
            film->addDebugLine(
                DebugLine(Vector2(sScreen.x, sScreen.y),
                Vector2(eScreen.x, eScreen.y)), debugRays[i].second);
        }
        const vector<pair<Vector3, Color> >& debugPoints = 
            debugData.getPoints();
        for(size_t i = 0; i < debugPoints.size(); ++i) {
            Vector3 pScreen = camera->worldToScreen(debugPoints[i].first);
            film->addDebugPoint(Vector2(pScreen.x, pScreen.y), 
                debugPoints[i].second);
        }        
    }
}
