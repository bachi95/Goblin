#include "GoblinBDPT.h"
#include "GoblinCamera.h"
#include "GoblinFilm.h"
#include "GoblinRay.h"

namespace Goblin {

    class BDPTTask : public RenderTask {
    public:
        BDPTTask(ImageTile* tile, BDPT* mLightTracer,
            const CameraPtr& camera, const ScenePtr& scene,
            const SampleQuota& sampleQuota, int samplePerPixel,
            int maxPathLength, RenderProgress* renderProgress);

        ~BDPTTask();

        void run();
    private:
        const BDPT* mBDPT;
        std::vector<PathVertex> mLightPath;
        std::vector<PathVertex> mEyePath;
        std::vector<BDPTMISNode> mMISNodes;
    };

    BDPTTask::BDPTTask(ImageTile* tile, BDPT* lightTracer,
        const CameraPtr& camera, const ScenePtr& scene,
        const SampleQuota& sampleQuota, int samplePerPixel,
        int maxPathLength, RenderProgress* renderProgress):
        RenderTask(tile, lightTracer, camera, scene, sampleQuota,
        samplePerPixel, renderProgress),
        mBDPT(lightTracer),
        mLightPath(maxPathLength + 1),
        mEyePath(maxPathLength + 1),
        mMISNodes(maxPathLength + 1) {}

    BDPTTask::~BDPTTask() {}

    void BDPTTask::run() {
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
                mBDPT->evalContribution(mScene, samples[s], *mRNG,
                    mLightPath, mEyePath, mMISNodes, mTile);
            }
            totalSampleCount += sampleNum;
        }
        mTile->setTotalSampleCount(totalSampleCount);
        delete [] samples;
        mRenderProgress->update();
    }

    BDPT::BDPT(int samplePerPixel, int threadNum,
        int maxPathLength, int debugS, int debugT, bool debugNoMIS):
        Renderer(samplePerPixel, threadNum),
        mTotalSamplesNum(0), mMaxPathLength(maxPathLength),
        mLightPathSampleIndexes(NULL), mEyePathSampleIndexes(NULL),
        mDebugS(debugS), mDebugT(debugT), mDebugNoMIS(debugNoMIS)
        {}

    BDPT::~BDPT() {
        if (mLightPathSampleIndexes) {
            delete [] mLightPathSampleIndexes;
            mLightPathSampleIndexes = NULL;
        }
        if (mEyePathSampleIndexes) {
            delete [] mEyePathSampleIndexes;
            mEyePathSampleIndexes = NULL;
        }
    }

    Color BDPT::Li(const ScenePtr& scene, const RayDifferential& ray,
        const Sample& sample, const RNG& rng,
        WorldDebugData* debugData) const {
        // bidirection path tracing doesn't use this camera based method
        return Color::Black;
    }

    void BDPT::evalContribution(const ScenePtr& scene,
        const Sample& sample, const RNG& rng,
        std::vector<PathVertex>& lightPath,
        std::vector<PathVertex>& eyePath,
        std::vector<BDPTMISNode>& misNodes,
        ImageTile* tile) const {
        // construct path randomwalk from light
        const vector<Light*>& lights = scene->getLights();
        if (lights.size() == 0) {
            return;
        }
        int lightVertexCount = constructLightPath(scene, sample, rng,
            lights, lightPath);
        // construct path randomwalk from camera
        const CameraPtr camera = scene->getCamera();
        int eyeVertexCount = constructEyePath(scene, sample, rng,
            camera, eyePath);
        for (int pathLength = 1; pathLength <= mMaxPathLength; ++pathLength) {
            int pathVertexCount = pathLength + 1;
            for (int s = 0; s <= pathVertexCount; ++s) {
                int t = pathVertexCount - s;
                // for debug purpose, only take certain strategy into account
                if (((mDebugS != -1) && (s != mDebugS)) ||
                    ((mDebugT != -1) && (t != mDebugT))) {
                    continue;
                }
                // we don't have long enough light or eye path for
                // this s/t combination case
                if (s > lightVertexCount || t > eyeVertexCount) {
                    continue;
                }
                // can not form a path with 0 or 1 vertex
                if ((s == 0 && t < 2) ||(t == 0 && s < 2) || (s + t) < 2) {
                    continue;
                }
                if ((t == 0 && !lightPath[s - 1].isCameraLens) ||
                    (s == 0 && !eyePath[t - 1].isLight())) {
                    continue;
                }
                // need to re-evaluate the contributed pixel coordinate for
                // t = 0 (light particle hit the camera lens) and
                // t = 1 (light path end vertex connect with camera lens) case
                Vector3 filmPixel(sample.imageX, sample.imageY, 0.0f);
                if (t == 0) {
                    filmPixel = camera->worldToScreen(
                        lightPath[s - 2].getPosition(),
                        lightPath[s - 1].getPosition());
                } else if (t == 1) {
                    filmPixel = camera->worldToScreen(
                        lightPath[s - 1].getPosition(),
                        eyePath[0].getPosition());
                }
                if (filmPixel == Camera::sInvalidPixel) {
                    continue;
                }
                // geometry factor between the light path end vertex and 
                // eye path end vertex
                float Gconnect = 1.0f;
                Color unweightedContribution = evalUnweightedContribution(
                    scene, camera, lightPath, s, eyePath, t, Gconnect);
                // no need to do the MIS evaluation if the path combination
                // has no contribution
                if (unweightedContribution == Color::Black) {
                    continue;
                }
                float weight = mDebugNoMIS ?
                    1.0f :
                    evalMIS(scene, camera, lightPath, s, eyePath, t,
                    Gconnect, misNodes);
                tile->addSample(filmPixel.x, filmPixel.y,
                    weight * unweightedContribution);
            }
        }
    }

    int BDPT::constructLightPath(const ScenePtr& scene, const Sample& sample,
        const RNG& rng, const std::vector<Light*>& lights,
        std::vector<PathVertex>& lightPath) const {
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
        float pdfBackward = pdfLightArea * pickLightPdf;
        float pdfLightDirection;
        BSDFSample bs(sample, mLightPathSampleIndexes[0], 0);
        Vector3 dir = light->sampleDirection(
            nLight, bs.uDirection[0], bs.uDirection[1], &pdfLightDirection);
        float pdfForward = pdfLightDirection / absdot(nLight, dir);
        lightPath[0] = PathVertex(Color(1.0f / pdfBackward),
            pLight, nLight, light, pdfForward, pdfBackward);
        Color throughput = lightPath[0].throughput / pdfForward;
        Ray ray(pLight, dir, 1e-3f);
        int lightVertexCount = 1;
        while (lightVertexCount <= mMaxPathLength) {
            float epsilon;
            Intersection isect;
            if (!scene->intersect(ray, &epsilon, &isect)) {
                break;
            }
            const Fragment& frag = isect.fragment;
            // for some light types (delta light like point/spot light),
            // we can't evaluate radiance before we know the intersection of
            // light particle shoot out from light
            if (lightVertexCount == 1) {
                throughput *= light->evalL(pLight, nLight, frag.getPosition());
            }
            BSDFSample bs(sample, mLightPathSampleIndexes[lightVertexCount], 0);
            Vector3 wo = -normalize(ray.d);
            Vector3 wi;
            float pdfW;
            BSDFType sampledType;
            Color f = isect.getMaterial()->sampleBSDF(frag, wo, bs,
                &wi, &pdfW, BSDFAll, &sampledType, BSDFImportance);
            bool isSpecular = (sampledType & BSDFSpecular) == BSDFSpecular;
            pdfForward = pdfW / absdot(wi, frag.getNormal());
            if (isSpecular) {
                pdfBackward = pdfForward;
            } else {
                pdfBackward = isect.getMaterial()->pdf(frag, wi, wo) /
                    absdot(wo, frag.getNormal());
            }
            lightPath[lightVertexCount] = PathVertex(throughput, isect,
                pdfForward, pdfBackward, isSpecular);
            lightPath[lightVertexCount].G = evalG(lightPath[lightVertexCount],
                lightPath[lightVertexCount - 1]);
            lightVertexCount += 1;
            if (f == Color::Black || pdfW == 0.0f) {
                break;
            }
            throughput *= f / pdfForward;
            ray = Ray(frag.getPosition(), wi, epsilon);
        }
        return lightVertexCount;
    }

    int BDPT::constructEyePath(const ScenePtr& scene,
        const Sample& sample, const RNG& rng,
        const CameraPtr& camera,
        std::vector<PathVertex>& eyePath) const {
        Vector3 nCamera;
        float pdfCamera;
        Vector3 pCamera = camera->samplePosition(sample, &nCamera, &pdfCamera);
        float pdfBackward = pdfCamera;
        float pdfEyeDirection;
        float We;
        Vector3 dir = camera->sampleDirection(
            sample, pCamera, &We, &pdfEyeDirection);
        float pdfForward = pdfEyeDirection / absdot(nCamera, dir);
        eyePath[0] = PathVertex(Color(1.0f / pdfBackward),
            pCamera, nCamera, camera.get(), pdfForward, pdfBackward);
        Color throughput = eyePath[0].throughput * We / pdfForward;
        Ray ray(pCamera, dir, 1e-3f);
        int eyeVertexCount = 1;
        while (eyeVertexCount <= mMaxPathLength) {
            float epsilon;
            Intersection isect;
            if (!scene->intersect(ray, &epsilon, &isect)) {
                break;
            }
            const Fragment& frag = isect.fragment;
            BSDFSample bs(sample, mEyePathSampleIndexes[eyeVertexCount], 0);
            Vector3 wo = -normalize(ray.d);
            Vector3 wi;
            float pdfW;
            BSDFType sampledType;
            Color f = isect.getMaterial()->sampleBSDF(frag, wo, bs,
                &wi, &pdfW, BSDFAll, &sampledType, BSDFRadiance);
            bool isSpecular = (sampledType & BSDFSpecular) == BSDFSpecular;
            pdfForward = pdfW / absdot(wi, frag.getNormal());
            if (isSpecular) {
                pdfBackward = pdfForward;
            } else {
                pdfBackward = isect.getMaterial()->pdf(frag, wi, wo) /
                    absdot(wo, frag.getNormal());
            }
            eyePath[eyeVertexCount] = PathVertex(throughput, isect,
                pdfForward, pdfBackward, isSpecular);
            eyePath[eyeVertexCount].G = evalG(eyePath[eyeVertexCount],
                eyePath[eyeVertexCount - 1]);
            eyeVertexCount += 1;
            if (f == Color::Black || pdfW == 0.0f) {
                break;
            }
            throughput *= f / pdfForward;
            ray = Ray(frag.getPosition(), wi, epsilon);
        }
        return eyeVertexCount;
    }

    Color BDPT::evalUnweightedContribution(
        const ScenePtr& scene, const CameraPtr& camera,
        const std::vector<PathVertex>& lightPath, int s,
        const std::vector<PathVertex>& eyePath, int t,
        float& G) const {
        // eval unweighted contribution
        Color aL = s == 0 ?
            Color::White : lightPath[s - 1].throughput;
        Color aE = t == 0 ?
            Color::White : eyePath[t - 1].throughput;
        Color cst;
        // eval connection factor (light path end point to eye path end point)
        if (s == 0) {
            cst = eyePath[t - 1].getLight()->evalL(
                eyePath[t - 1].getPosition(),
                eyePath[t - 1].getNormal(),
                eyePath[t - 2].getPosition());
        } else if (t == 0) {
            cst = Color(camera->evalWe(lightPath[s - 1].getPosition(),
                lightPath[s - 2].getPosition()));
        } else {
            const PathVertex& sEndV = lightPath[s - 1];
            const PathVertex& tEndV = eyePath[t - 1]; 
            Vector3 connectVector = tEndV.getPosition() - sEndV.getPosition();
            Vector3 connectDir = normalize(connectVector);
            Color fsL;
            if (s == 1) {
                fsL = sEndV.light->evalL(sEndV.getPosition(),
                    sEndV.getNormal(), tEndV.getPosition());
            } else {
                const Vector3 wi = connectDir;
                const Vector3 wo = normalize(lightPath[s - 2].getPosition() -
                    sEndV.getPosition());
                fsL = lightPath[s - 1].material->bsdf(sEndV.fragment, wo, wi,
                    BSDFAll, BSDFImportance);
            }
            if (fsL == Color::Black) {
                return Color::Black;
            }

            Color fsE;
            if (t == 1) {
                fsE = Color(camera->evalWe(tEndV.getPosition(),
                    sEndV.getPosition()));
            } else {
                const Vector3 wi = -connectDir;
                const Vector3 wo = normalize(
                    eyePath[t - 2].getPosition() -
                    tEndV.getPosition());
                fsE = tEndV.material->bsdf(tEndV.fragment, wo, wi,
                    BSDFAll, BSDFRadiance);
            }
            if (fsE == Color::Black) {
                return Color::Black;
            }
            G = evalG(sEndV, tEndV);
            if (G == 0.0f){
                return Color::Black;
            }
            // occlude test
            float occludeDistance = length(connectVector);
            float epsilon = 1e-3f * occludeDistance;
            Ray occludeRay(sEndV.getPosition(), connectDir,
                epsilon, occludeDistance - epsilon);
            if (scene->intersect(occludeRay)) {
                return Color::Black;
            }
            cst = fsL * G * fsE;
        }
        return aL * cst * aE;
    }

    float BDPT::evalG(const PathVertex& a, const PathVertex& b) const {
        const Vector3& pA = a.getPosition();
        const Vector3& pB = b.getPosition();
        const Vector3 vAB = pB - pA;
        float invLengthAB = 1.0f / length(vAB);
        const Vector3 dAB = vAB * invLengthAB;
        // TODO right now all the camera types return a valid normal
        // (eye face toward film) so we don't treat delta camera as
        // exceptional case. We need to handle further exceptional case
        // when adding support of IBL
        // (both invLength and cos need special treatment)
        bool deltaExceptionA = (a.isLight() && a.getLight()->isDelta());
        float cosA = deltaExceptionA ? 1.0f : absdot(a.getNormal(), dAB);
        bool deltaExceptionB = (b.isLight() && b.getLight()->isDelta());
        float cosB = deltaExceptionB ? 1.0f : absdot(b.getNormal(), dAB);
        return cosA * cosB * invLengthAB * invLengthAB;
    }

    float BDPT::evalMIS(const ScenePtr& scene, const CameraPtr& camera,
        const std::vector<PathVertex>& lightPath, int s,
        const std::vector<PathVertex>& eyePath, int t,
        const float GConnect, std::vector<BDPTMISNode>& misNodes) const {
        // re-evaluate the pdfForward, pdfBackward for lightPath and eyePath
        // end vertex
        float pdfSEndForward, pdfSEndBackward;
        float pdfTEndForward, pdfTEndBackward;
        if (s == 0) {
            // eye path end vertex is a light
            const Vector3& p = eyePath[t - 1].getPosition();
            const Vector3& n = eyePath[t - 1].getNormal();
            const Light* light = eyePath[t - 1].light;
            pdfTEndForward = mPickLightPdf[light->getId()] *
                light->pdfPosition(scene, p);
            const Vector3& wo = normalize(eyePath[t - 2].getPosition() - p);
            pdfTEndBackward = light->pdfDirection(p, n, wo) / dot(n, wo);
            // this two are not used in s = 0 case. Set to 0 to be explicit.
            pdfSEndForward = 0.0f;
            pdfSEndBackward = 0.0f;
        } else if (t == 0) {
            // light path end vertex is camera lens
            const Vector3& p = lightPath[s - 1].getPosition();
            const Vector3& n = lightPath[s - 1].getNormal();
            pdfSEndForward = camera->pdfPosition(p);
            const Vector3 wo = normalize(lightPath[s - 2].getPosition() - p);
            pdfSEndBackward = camera->pdfDirection(p, wo) / dot(n, wo);
            // this two are not used in t == 0 case. Set to 0 to be explicit.
            pdfTEndForward = 0.0f;
            pdfTEndBackward = 0.0f;
        } else {
            const PathVertex& sEnd = lightPath[s - 1];
            const PathVertex& tEnd = eyePath[t - 1];
            const Vector3& pSEnd = sEnd.getPosition();
            const Vector3& nSEnd = sEnd.getNormal();
            const Vector3& pTEnd = tEnd.getPosition();
            const Vector3& nTEnd = tEnd.getNormal();
            const Vector3 dSToT = normalize(pTEnd - pSEnd);
            if (s == 1) {
                float pdfW = sEnd.light->pdfDirection(pSEnd, nSEnd, dSToT);
                pdfSEndForward = pdfW / dot(nSEnd, dSToT);
                pdfSEndBackward = sEnd.pdfBackward;
            } else {
                const Vector3 dSEndToSPrev = normalize(
                    lightPath[s - 2].getPosition() - pSEnd);
                pdfSEndForward = sEnd.material->pdf(sEnd.fragment,
                    dSEndToSPrev, dSToT) / dot(dSToT, nSEnd);
                pdfSEndBackward = sEnd.material->pdf(sEnd.fragment,
                    dSToT, dSEndToSPrev) / dot(dSEndToSPrev, nSEnd);
            }
            const Vector3 dTToS = -dSToT;
            if (t == 1) {
                float pdfW = camera->pdfDirection(pTEnd, dTToS);
                pdfTEndForward =  pdfW / dot(nTEnd, dTToS);
                pdfTEndBackward = tEnd.pdfBackward;
            } else {
                const Vector3 dTEndToTPrev = normalize(
                    eyePath[t - 2].getPosition() - pTEnd);
                pdfTEndForward = tEnd.material->pdf(tEnd.fragment,
                    dTEndToTPrev, dTToS) / dot(dTToS, nTEnd);
                pdfTEndBackward = tEnd.material->pdf(tEnd.fragment,
                    dTToS, dTEndToTPrev) / dot(dTEndToTPrev, nTEnd);
            }
        }

        // fill the pdf in lightPath and eyePath to misNodes
        int k = s + t - 1;
        for (int i = 0; i < s - 1; ++i) {
            misNodes[i].pTowardLight = i == 0 ?
                lightPath[0].pdfBackward :
                lightPath[i].pdfBackward * lightPath[i].G;
            misNodes[i].pTowardEye =
                lightPath[i].pdfForward * lightPath[i + 1].G;
            misNodes[i].isSpecular = lightPath[i].isSpecular;
        }
        if (s > 0) {
            misNodes[s - 1].pTowardLight = s == 1 ?
                pdfSEndBackward : pdfSEndBackward * lightPath[s - 1].G;
            misNodes[s - 1].pTowardEye = ((s - 1) == k) ?
                pdfSEndForward : pdfSEndForward * GConnect;
            misNodes[s - 1].isSpecular = lightPath[s - 1].isSpecular;
        }
        for (int i = 0; i < t - 1; ++i) {
            misNodes[k - i].pTowardEye = i == 0 ?
                eyePath[0].pdfBackward :
                eyePath[i].pdfBackward * eyePath[i].G;
            misNodes[k - i].pTowardLight =
                eyePath[i].pdfForward * eyePath[i + 1].G;
            misNodes[k - i].isSpecular = eyePath[i].isSpecular;
        }
        if (t > 0) {
            misNodes[k - (t - 1)].pTowardEye = t == 1 ?
                pdfTEndBackward : pdfTEndBackward * eyePath[t - 1].G;
            misNodes[k - (t - 1)].pTowardLight = ((t - 1) == k) ?
                pdfTEndForward : pdfTEndForward * GConnect;
            misNodes[k - (t - 1)].isSpecular = eyePath[t - 1].isSpecular;
        }

        // iterate from the connection end point to calculate
        // relative pdfA and add it to misWeightSum (power heuristic)
        float pK = 1.0f;
        float misWeightSum = 1.0f;
        for (int i = s; i <= k; ++i) {
            if (i == 0) {
                pK *= misNodes[0].pTowardLight / misNodes[1].pTowardLight;
                // exception handling for specular case
                if (misNodes[1].isSpecular) {
                    continue;
                }
            } else if (i == k) {
                // there is no way that light path can hit a delta camera
                if (camera->isDelta()) {
                    break;
                }
                pK *= misNodes[k - 1].pTowardEye / misNodes[k].pTowardEye;
            } else {
                pK *= misNodes[i - 1].pTowardEye / misNodes[i + 1].pTowardLight;
                // exception handling for specular case
                if (misNodes[i].isSpecular || misNodes[i + 1].isSpecular) {
                    continue;
                }
            }
            misWeightSum += pK * pK;
        }
        pK = 1.0f;
        for (int i = s; i > 0; --i) {
            if (i == (k + 1)) {
                pK *= misNodes[k].pTowardEye / misNodes[k - 1].pTowardEye;
                // exception handling for specular case
                if (misNodes[k - 1].isSpecular) {
                    continue;
                }
            } else if (i == 1) {
                // thre is no way that eye path can hit a delta light
                if (lightPath[0].light->isDelta()) {
                    break;
                }
                pK *= misNodes[1].pTowardLight / misNodes[0].pTowardLight;
            } else {
                pK *= misNodes[i].pTowardLight / misNodes[i - 2].pTowardEye;
                // exception handling for specular case
                if (misNodes[i - 1].isSpecular || misNodes[i - 2].isSpecular) {
                    continue;
                }
            }
            misWeightSum += pK * pK;
        }

        return 1.0f / misWeightSum;
    }

    void BDPT::render(const ScenePtr& scene) {
        const CameraPtr camera = scene->getCamera();
        Film* film = camera->getFilm();
        SampleQuota sampleQuota;
        querySampleQuota(scene, &sampleQuota);

        vector<ImageTile*>& tiles = film->getTiles();
        vector<Task*> bdptTasks;
        RenderProgress progress((int)tiles.size());
        for (size_t i = 0; i < tiles.size(); ++i) {
            bdptTasks.push_back(new BDPTTask(tiles[i], this,
                camera, scene, sampleQuota, mSamplePerPixel,
                mMaxPathLength, &progress));
        }

        ThreadPool threadPool(mThreadNum);
        threadPool.enqueue(bdptTasks);
        threadPool.waitForAll();
        //clean up
        for (size_t i = 0; i < bdptTasks.size(); ++i) {
            delete bdptTasks[i];
        }
        bdptTasks.clear();

        uint64_t totalSampleCount = 0;
        for (size_t i = 0; i < tiles.size(); ++i) {
            totalSampleCount += tiles[i]->getTotalSampleCount();
        }
        film->mergeTiles();
        film->scaleImage(film->getFilmArea() / totalSampleCount);
        film->writeImage(false);
    }

    void BDPT::querySampleQuota(const ScenePtr& scene,
        SampleQuota* sampleQuota){

        if (mLightSampleIndexes) {
            delete [] mLightSampleIndexes;
            mLightSampleIndexes = NULL;
        }
        if (mLightPathSampleIndexes) {
            delete [] mLightPathSampleIndexes;
            mLightPathSampleIndexes = NULL;
        }
        if (mEyePathSampleIndexes) {
            delete [] mEyePathSampleIndexes;
            mEyePathSampleIndexes = NULL;
        }
        if (mPowerDistribution) {
            delete mPowerDistribution;
            mPowerDistribution = NULL;
        }
        if (mPickLightSampleIndexes) {
            delete [] mPickLightSampleIndexes;
            mPickLightSampleIndexes = NULL;
        }

        mLightSampleIndexes = new LightSampleIndex[1];
        mLightSampleIndexes[0] = LightSampleIndex(sampleQuota, 1);
        mPickLightSampleIndexes = new SampleIndex[1];
        mPickLightSampleIndexes[0] = sampleQuota->requestOneDQuota(1);

        mLightPathSampleIndexes = new BSDFSampleIndex[mMaxPathLength + 1];
        for (int i = 0; i < mMaxPathLength + 1; ++i) {
            mLightPathSampleIndexes[i] = BSDFSampleIndex(sampleQuota, 1);
        }

        mEyePathSampleIndexes = new BSDFSampleIndex[mMaxPathLength + 1];
        for (int i = 0; i < mMaxPathLength + 1; ++i) {
            mEyePathSampleIndexes[i] = BSDFSampleIndex(sampleQuota, 1);
        }

        const vector<Light*>& lights = scene->getLights();
        std::vector<float> lightPowers;
        size_t maxLightId = 0;
        float totalPower = 0.0f;
        for (size_t i = 0; i < lights.size(); ++i) {
            float power = lights[i]->power(scene).luminance();
            lightPowers.push_back(power);
            totalPower += power;
            if (lights[i]->getId() > maxLightId) {
                maxLightId = lights[i]->getId();
            }
        }
        mPowerDistribution = new CDF1D(lightPowers);
        // build a id -> pick light pdf table since we need to evaluate
        // the pdf to pick a specified light when we are doing the MIS
        // evaluation
        mPickLightPdf.resize(maxLightId + 1, 0.0f);
        for (size_t i = 0; i < lights.size(); ++i) {
            mPickLightPdf[lights[i]->getId()] = (lightPowers[i] / totalPower);
        }
    }

    Renderer* BDPTCreator::create(const ParamSet& params) const {
        int samplePerPixel = params.getInt("sample_per_pixel", 1);
        int threadNum = params.getInt("thread_num",
            boost::thread::hardware_concurrency());
        int maxPathLength = params.getInt("max_path_length", 5);
        int debugS = params.getInt("debug_s", -1);
        int debugT = params.getInt("debug_t", -1);
        bool debugNoMIS = params.getBool("debug_no_mis", false);
        return new BDPT(samplePerPixel, threadNum, maxPathLength,
            debugS, debugT, debugNoMIS);
    }
}
