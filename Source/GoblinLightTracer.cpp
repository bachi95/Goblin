#include "GoblinLightTracer.h"
#include "GoblinCamera.h"
#include "GoblinFilm.h"
#include "GoblinRay.h"

namespace Goblin {

    LightTracer::LightTracer(int samplePerPixel, int maxPathLength):
        Renderer(samplePerPixel, 1),
        mTotalSamplesNum(0), mMaxPathLength(maxPathLength),
        mPathVertices(mMaxPathLength){}

    LightTracer::~LightTracer() {}

    Color LightTracer::Li(const ScenePtr& scene, const RayDifferential& ray,
        const Sample& sample, const RNG& rng,
        WorldDebugData* debugData) const {
        // light tracing method don't use this camera based method actually
        return Color::Black;
    }

    void LightTracer::splatFilm(const ScenePtr& scene, const Sample& sample,
        const RNG& rng) {
        const vector<Light*>& lights = scene->getLights();
        if (lights.size() == 0) {
            return;
        }
        // get the camera point
        const CameraPtr camera = scene->getCamera();
        Film* film = camera->getFilm();
        Vector3 nCamera;
        float pdfCamera;
        Vector3 pCamera = camera->samplePosition(sample, &nCamera, &pdfCamera);
        PathVertex cVertex(Color(1.0f / pdfCamera), pCamera,
            nCamera, camera.get());

        // get the light point
        float pickLightPdf;
        float pickSample = sample.u1D[mPickLightSampleIndexes[0].offset][0];
        int lightIndex = mPowerDistribution->sampleDiscrete(
            pickSample, &pickLightPdf);
        const Light* light = lights[lightIndex];
        LightSample ls(sample, mLightSampleIndexes[0], 0);
        Vector3 nLight;
        float pdfLightArea;
        Vector3 pLight = light->samplePosition(scene, ls, &nLight,
            &pdfLightArea);

        mPathVertices[0] = PathVertex(
            Color(1.0f / (pdfLightArea * pickLightPdf)),
            pLight, nLight, light);
        float pdfLightDirection;
        BSDFSample bs(sample, mBSDFSampleIndexes[0], 0);
        Vector3 dir = light->sampleDirection(
            nLight, bs.uDirection[0], bs.uDirection[1], &pdfLightDirection);
        Color throughput = mPathVertices[0].throughput *
            absdot(nLight, dir) / pdfLightDirection;
        Ray ray(pLight, dir, 1e-5f);
        int lightVertex = 1;
        for (; lightVertex < mMaxPathLength; ++lightVertex) {
            float epsilon;
            Intersection isect;
            if (!scene->intersect(ray, &epsilon, &isect)) {
                break;
            }
            const Fragment& frag = isect.fragment;
            mPathVertices[lightVertex] = PathVertex(throughput, frag,
                isect.getMaterial().get());
            BSDFSample bs(sample, mBSDFSampleIndexes[lightVertex], 0);
            Vector3 wo = -normalize(ray.d);
            Vector3 wi;
            float pdfW;
            // we should add an adjoint option to bsdf, for now
            // we are testing only on lambert case so it should
            // be fine....
            Color f = isect.getMaterial()->sampleBSDF(frag, wo, bs,
                &wi, &pdfW);
            throughput *= f * absdot(wi, frag.getNormal()) / pdfW;
            ray.o = frag.getPosition();
            ray.d = wi;
            ray.mint = 1e-4f;
        }

        for (int s = 1; s <= lightVertex; ++s) {

            const PathVertex& pv = mPathVertices[s  - 1];
            const Vector3& pvPos = pv.fragment.getPosition();
            Vector3 filmPixel = camera->worldToScreen(
                pv.fragment.getPosition(), pCamera);
            if (filmPixel == Camera::sInvalidPixel) {
                continue;
            }
            // occlude test
            Vector3 pv2Cam(pCamera - pvPos);
            Ray occludeRay(pvPos, normalize(pv2Cam),
                1e-5f, length(pv2Cam));
            if (scene->intersect(occludeRay)) {
                continue;
            }
            Color fsL;
            Vector3 wo = normalize(pCamera - pvPos);
            if (s > 1) {
                Vector3 wi =
                    normalize(mPathVertices[s - 2].fragment.getPosition() - pvPos);
                Color f = pv.material->bsdf(pv.fragment, wo, wi);
                fsL = f * light->evalL(pLight, nLight,
                    mPathVertices[1].fragment.getPosition());
            } else {
                fsL = pv.light->evalL(pLight, nLight, pCamera);
            }
            float fsE = camera->evalWe(pCamera, pvPos);
            float G = absdot(pv.fragment.getNormal(), wo) * absdot(nCamera, wo) /
                squaredLength(pvPos - pCamera);
            
            Color pathContribution = fsL * fsE * G *
                pv.throughput * cVertex.throughput;
            film->addSample(filmPixel.x, filmPixel.y, pathContribution);
        }
    }

    void LightTracer::render(const ScenePtr& scene) {
        const CameraPtr camera = scene->getCamera();
        Film* film = camera->getFilm();
        SampleQuota sampleQuota;
        querySampleQuota(scene, &sampleQuota);
        int xStart, xEnd, yStart, yEnd;
        film->getSampleRange(&xStart, &xEnd, &yStart, &yEnd);
        RNG rng;
        Sampler sampler(xStart, xEnd, yStart, yEnd, mSamplePerPixel,
            sampleQuota, &rng);
        int batchAmount = sampler.maxSamplesPerRequest();
        Sample* samples = sampler.allocateSampleBuffer(batchAmount);
        mTotalSamplesNum = sampler.maxTotalSamples();
        cout <<"total samples N: " <<mTotalSamplesNum << endl;
        int sampleNum = 0;
        while((sampleNum = sampler.requestSamples(samples)) > 0) {
            for (int s = 0; s <sampleNum; ++s) {
                splatFilm(scene, samples[s], rng);
            }
        }
        delete [] samples;
        film->scaleImage(film->getFilmArea() / mTotalSamplesNum);
        film->writeImage(false);
    }

    void LightTracer::querySampleQuota(const ScenePtr& scene,
        SampleQuota* sampleQuota){

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

        mLightSampleIndexes = new LightSampleIndex[1];
        mLightSampleIndexes[0] = LightSampleIndex(sampleQuota, 1);
        mPickLightSampleIndexes = new SampleIndex[1];
        mPickLightSampleIndexes[0] = sampleQuota->requestOneDQuota(1);
        mBSDFSampleIndexes = new BSDFSampleIndex[mMaxPathLength];
        for(int i = 0; i < mMaxPathLength; ++i) {
            mBSDFSampleIndexes[i] = BSDFSampleIndex(sampleQuota, 1);
        }

        const vector<Light*>& lights = scene->getLights();
        vector<float> lightPowers;
        for(size_t i = 0; i < lights.size(); ++i) {
            lightPowers.push_back(lights[i]->power(scene).luminance());
        }
        mPowerDistribution = new CDF1D(lightPowers);
    }

    Renderer* LightTracerCreator::create(const ParamSet& params) const {
        int samplePerPixel = params.getInt("sample_per_pixel", 1);
        int maxPathLength = params.getInt("max_path_length", 5);
        return new LightTracer(samplePerPixel, maxPathLength);
    }
}
