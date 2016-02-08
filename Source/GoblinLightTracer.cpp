#include "GoblinLightTracer.h"
#include "GoblinCamera.h"
#include "GoblinFilm.h"
#include "GoblinRay.h"

namespace Goblin {

    class LightTraceTask : public RenderTask {
    public:
        LightTraceTask(ImageTile* tile, LightTracer* mLightTracer,
            const CameraPtr& camera, const ScenePtr& scene,
            const SampleQuota& sampleQuota, int samplePerPixel,
            int maxPathLength, RenderProgress* renderProgress);
        ~LightTraceTask();
        void run();
    private:
        const LightTracer* mLightTracer;
        std::vector<PathVertex> mPathVertices;
    };

    LightTraceTask::LightTraceTask(ImageTile* tile, LightTracer* lightTracer,
        const CameraPtr& camera, const ScenePtr& scene,
        const SampleQuota& sampleQuota, int samplePerPixel,
        int maxPathLength, RenderProgress* renderProgress):
        RenderTask(tile, lightTracer, camera, scene, sampleQuota,
        samplePerPixel, renderProgress),
        mLightTracer(lightTracer), mPathVertices(maxPathLength + 1) {}

    LightTraceTask::~LightTraceTask() {}

    void LightTraceTask::run() {
        int xStart, xEnd, yStart, yEnd;
        mTile->getSampleRange(&xStart, &xEnd, &yStart, &yEnd);
        Sampler sampler(xStart, xEnd, yStart, yEnd,
            mSamplePerPixel, mSampleQuota, mRNG);
        int batchAmount = sampler.maxSamplesPerRequest();
        Sample* samples = sampler.allocateSampleBuffer(batchAmount);
        WorldDebugData debugData;
        int sampleNum = 0;
        uint64_t totalSampleCount = 0;
        while ((sampleNum = sampler.requestSamples(samples)) > 0) {
            for (int s = 0; s <sampleNum; ++s) {
                //mLightTracer->splatFilmT0(mScene, samples[s], *mRNG,
                //    mPathVertices, mTile);

                mLightTracer->splatFilmT1(mScene, samples[s], *mRNG,
                    mPathVertices, mTile);

                //mLightTracer->splatFilmS1(mScene, samples[s], *mRNG,
                //    mPathVertices, mTile);
            }
            totalSampleCount += sampleNum;
        }
        mTile->setTotalSampleCount(totalSampleCount);
        delete [] samples;
        mRenderProgress->update();
    }

    LightTracer::LightTracer(int samplePerPixel, int threadNum,
        int maxPathLength):
        Renderer(samplePerPixel, threadNum),
        mTotalSamplesNum(0), mMaxPathLength(maxPathLength) {}

    LightTracer::~LightTracer() {}

    Color LightTracer::Li(const ScenePtr& scene, const RayDifferential& ray,
        const Sample& sample, const RNG& rng,
        WorldDebugData* debugData) const {
        // light tracing method doesn't use this camera based method actually
        return Color::Black;
    }

    void LightTracer::splatFilmT1(const ScenePtr& scene, const Sample& sample,
        const RNG& rng, std::vector<PathVertex>& pathVertices,
        ImageTile* tile) const {
        const vector<Light*>& lights = scene->getLights();
        if (lights.size() == 0) {
            return;
        }
        // get the camera point
        const CameraPtr camera = scene->getCamera();
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

        pathVertices[0] = PathVertex(
            Color(1.0f / (pdfLightArea * pickLightPdf)),
            pLight, nLight, light);
        float pdfLightDirection;
        BSDFSample bs(sample, mBSDFSampleIndexes[0], 0);
        Vector3 dir = light->sampleDirection(
            nLight, bs.uDirection[0], bs.uDirection[1], &pdfLightDirection);
        Color throughput = pathVertices[0].throughput *
            absdot(nLight, dir) / pdfLightDirection;
        Ray ray(pLight, dir, 1e-3f);
        int lightVertex = 1;
        while (lightVertex < mMaxPathLength) {
            float epsilon;
            Intersection isect;
            if (!scene->intersect(ray, &epsilon, &isect)) {
                break;
            }
            const Fragment& frag = isect.fragment;
            pathVertices[lightVertex] = PathVertex(throughput, isect);
            BSDFSample bs(sample, mBSDFSampleIndexes[lightVertex], 0);
            lightVertex += 1;
            Vector3 wo = -normalize(ray.d);
            Vector3 wi;
            float pdfW;
            Color f = isect.getMaterial()->sampleBSDF(frag, wo, bs,
                &wi, &pdfW, BSDFAll, NULL, BSDFImportance);
            if (f == Color::Black || pdfW == 0.0f) {
                break;
            }
            throughput *= f * absdot(wi, frag.getNormal()) / pdfW;
            ray = Ray(frag.getPosition(), wi, epsilon);
        }

        for (int s = 1; s <= lightVertex; ++s) {
            const PathVertex& pv = pathVertices[s  - 1];
            const Vector3& pvPos = pv.fragment.getPosition();
            Vector3 filmPixel = camera->worldToScreen(
                pv.fragment.getPosition(), pCamera);
            if (filmPixel == Camera::sInvalidPixel) {
                continue;
            }
            // occlude test
            Vector3 pv2Cam(pCamera - pvPos);
            float occludeDistance = length(pv2Cam);
            float epsilon = 1e-3f * occludeDistance;
            Ray occludeRay(pvPos, normalize(pv2Cam),
                epsilon, occludeDistance - epsilon);
            if (scene->intersect(occludeRay)) {
                continue;
            }
            Color fsL;
            Vector3 wo = normalize(pCamera - pvPos);
            if (s > 1) {
                Vector3 wi =
                    normalize(pathVertices[s - 2].fragment.getPosition() - pvPos);
                Color f = pv.material->bsdf(pv.fragment, wo, wi,
                    BSDFAll, BSDFImportance);
                fsL = f * light->evalL(pLight, nLight,
                    pathVertices[1].fragment.getPosition());
            } else {
                fsL = pv.light->evalL(pLight, nLight, pCamera);
            }
            float fsE = camera->evalWe(pCamera, pvPos);
            float G = absdot(pv.fragment.getNormal(), wo) * absdot(nCamera, wo) /
                squaredLength(pvPos - pCamera);
            
            Color pathContribution = fsL * fsE * G *
                pv.throughput * cVertex.throughput;
            tile->addSample(filmPixel.x, filmPixel.y, pathContribution);
        }
    }

    void LightTracer::splatFilmT0(const ScenePtr& scene, const Sample& sample,
        const RNG& rng, std::vector<PathVertex>& pathVertices,
        ImageTile* tile) const {
        const vector<Light*>& lights = scene->getLights();
        if (lights.size() == 0) {
            return;
        }
        const CameraPtr camera = scene->getCamera();
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

        pathVertices[0] = PathVertex(
            Color(1.0f / (pdfLightArea * pickLightPdf)),
            pLight, nLight, light);
        float pdfLightDirection;
        BSDFSample bs(sample, mBSDFSampleIndexes[0], 0);
        Vector3 dir = light->sampleDirection(
            nLight, bs.uDirection[0], bs.uDirection[1], &pdfLightDirection);
        Color throughput = pathVertices[0].throughput *
            absdot(nLight, dir) / pdfLightDirection;
        Ray ray(pLight, dir, 1e-3f);
        int lightVertex = 1;
        while (lightVertex <= mMaxPathLength) {
            float epsilon;
            Intersection isect;
            if (!scene->intersect(ray, &epsilon, &isect)) {
                break;
            }
            const Fragment& frag = isect.fragment;
            pathVertices[lightVertex] = PathVertex(throughput, isect);
            // last light vertex hit the camera lens, evaluate result and done
            if (isect.isCameraLens()) {
                const Vector3& pCamera = frag.getPosition();
                const Vector3& pS_1 =
                    pathVertices[lightVertex - 1].fragment.getPosition();
                Vector3 filmPixel = camera->worldToScreen(pS_1, pCamera);
                if (filmPixel != Camera::sInvalidPixel) {
                    Color L = light->evalL(pLight, nLight,
                        pathVertices[1].fragment.getPosition());
                    float We = camera->evalWe(pCamera, pS_1);
                    tile->addSample(filmPixel.x, filmPixel.y,
                        L * We * pathVertices[lightVertex].throughput);
                }
                break;
            }
            // continue random walk the light vertex path until it hits lens
            BSDFSample bs(sample, mBSDFSampleIndexes[lightVertex], 0);
            lightVertex += 1;
            Vector3 wo = -normalize(ray.d);
            Vector3 wi;
            float pdfW;
            Color f = isect.getMaterial()->sampleBSDF(frag, wo, bs,
                &wi, &pdfW, BSDFAll, NULL, BSDFImportance);
            if (f == Color::Black || pdfW == 0.0f) {
                break;
            }
            throughput *= f * absdot(wi, frag.getNormal()) / pdfW;
            ray = Ray(frag.getPosition(), wi, epsilon);
        }
    }

    void LightTracer::splatFilmS1(const ScenePtr& scene, const Sample& sample,
        const RNG& rng, std::vector<PathVertex>& pathVertices,
        ImageTile* tile) const {
        const vector<Light*>& lights = scene->getLights();
        if (lights.size() == 0) {
            return;
        }
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
        PathVertex lVertex(Color(1.0f / (pdfLightArea * pickLightPdf)),
            pLight, nLight, light);

        // get the camera point
        const CameraPtr camera = scene->getCamera();
        Vector3 nCamera;
        float pdfCamera;
        Vector3 pCamera = camera->samplePosition(sample, &nCamera, &pdfCamera);

        pathVertices[0] = PathVertex(Color(1.0f / pdfCamera),
            pCamera, nCamera, camera.get());

        float pdfEyeDirection;
        float We;
        Vector3 dir = camera->sampleDirection(
            sample, pCamera, &We, &pdfEyeDirection);
        Color throughput = pathVertices[0].throughput *
            absdot(nCamera, dir) / pdfEyeDirection;
        Ray ray(pCamera, dir, 1e-3f);
        int eyeVertex = 1;
        while (eyeVertex < mMaxPathLength) {
            float epsilon;
            Intersection isect;
            if (!scene->intersect(ray, &epsilon, &isect)) {
                break;
            }
            const Fragment& frag = isect.fragment;
            pathVertices[eyeVertex] = PathVertex(throughput, isect);
            BSDFSample bs(sample, mBSDFSampleIndexes[eyeVertex], 0);
            eyeVertex += 1;
            Vector3 wo = -normalize(ray.d);
            Vector3 wi;
            float pdfW;
            Color f = isect.getMaterial()->sampleBSDF(frag, wo, bs,
                &wi, &pdfW);
            if (f == Color::Black || pdfW == 0.0f) {
                break;
            }
            throughput *= f * absdot(wi, frag.getNormal()) / pdfW;
            ray = Ray(frag.getPosition(), wi, epsilon);
        }

        for (int t = 1; t <= eyeVertex; ++t) {
            const PathVertex& pv = pathVertices[t  - 1];
            const Vector3& pvPos = pv.fragment.getPosition();
            // occlude test
            Vector3 pv2Light(pLight - pvPos);
            float occludeDistance = length(pv2Light);
            float epsilon = 1e-3f * occludeDistance;
            Ray occludeRay(pvPos, normalize(pv2Light),
                epsilon, occludeDistance - epsilon);
            if (scene->intersect(occludeRay)) {
                continue;
            }
            Color fsE;
            Vector3 filmPixel;
            Vector3 wi = normalize(pLight - pvPos);
            if (t > 1) {
                Vector3 wo =
                    normalize(pathVertices[t - 2].fragment.getPosition() - pvPos);
                Color f = pv.material->bsdf(pv.fragment, wo, wi);
                fsE = f * We;
                filmPixel.x = sample.imageX;
                filmPixel.y = sample.imageY;
            } else {
                filmPixel = camera->worldToScreen(pLight, pvPos);
                fsE = Color(camera->evalWe(pCamera, pLight));
            }
            Color fsL = light->evalL(pLight, nLight, pvPos);
            float G = absdot(pv.fragment.getNormal(), wi) * absdot(nLight, wi) /
                squaredLength(pLight - pvPos);
            Color pathContribution = fsL * fsE * G *
                pv.throughput * lVertex.throughput;
            tile->addSample(filmPixel.x, filmPixel.y, pathContribution);
        }
    }

    void LightTracer::render(const ScenePtr& scene) {
        const CameraPtr camera = scene->getCamera();
        Film* film = camera->getFilm();
        SampleQuota sampleQuota;
        querySampleQuota(scene, &sampleQuota);

        vector<ImageTile*>& tiles = film->getTiles();
        vector<Task*> lightTraceTasks;
        RenderProgress progress((int)tiles.size());
        for(size_t i = 0; i < tiles.size(); ++i) {
            lightTraceTasks.push_back(new LightTraceTask(tiles[i], this,
                camera, scene, sampleQuota, mSamplePerPixel,
                mMaxPathLength, &progress));
        }

        ThreadPool threadPool(mThreadNum);
        threadPool.enqueue(lightTraceTasks);
        threadPool.waitForAll();
        //clean up
        for(size_t i = 0; i < lightTraceTasks.size(); ++i) {
            delete lightTraceTasks[i];
        }
        lightTraceTasks.clear();

        uint64_t totalSampleCount = 0;
        for (size_t i = 0; i < tiles.size(); ++i) {
            totalSampleCount += tiles[i]->getTotalSampleCount();
        }
        film->mergeTiles();
        film->scaleImage(film->getFilmArea() / totalSampleCount);
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
        mBSDFSampleIndexes = new BSDFSampleIndex[mMaxPathLength + 1];
        for(int i = 0; i < mMaxPathLength + 1; ++i) {
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
        int threadNum = params.getInt("thread_num",
            boost::thread::hardware_concurrency());
        int maxPathLength = params.getInt("max_path_length", 5);
        return new LightTracer(samplePerPixel, threadNum, maxPathLength);
    }
}
