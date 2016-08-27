#include "GoblinSPPM.h"
#include "GoblinCamera.h"
#include "GoblinFilm.h"
#include "GoblinRay.h"

namespace Goblin {

    const float PixelData::sInvalidRadius = -1.0f;

    class RayTraceTLS : public ThreadLocalStorage {
    public:
        RayTraceTLS(const SampleQuota& sampleQuota) {
            mSample.allocateQuota(sampleQuota);
        }
        Sample mSample;
    };

    class RayTraceTLSManager : public TLSManager {
    public:
        RayTraceTLSManager(const SampleQuota& sampleQuota):
            mSampleQuota(sampleQuota) {}

        void initialize(TLSPtr& tlsPtr) {
            tlsPtr.reset(new RayTraceTLS(mSampleQuota));
        }

        void finalize(TLSPtr& tlsPtr) {}
    private:
        const SampleQuota& mSampleQuota;
    };


    class RayTraceTask : public Task {
    public:
        RayTraceTask(SPPM* sppm, const ScenePtr& scene,
            const SampleRange& sampleRange,
            const PermutedHalton& halton):
            mSPPM(sppm), mScene(scene), mCurrentIteration(0),
            mSampleRange(sampleRange), mHalton(halton) {
            // make each pixel uses different QMC sub sequence
            mHaltonStartID.resize((mSampleRange.xEnd - mSampleRange.xStart) *
                (mSampleRange.yEnd - mSampleRange.yStart));
            for (size_t i = 0; i < mHaltonStartID.size(); ++i) {
                mHaltonStartID[i] = mRNG.randomUInt();
            }
        }

        void run(TLSPtr& tls);

        void nextIteration() { mCurrentIteration++; }
    private:
        SPPM* mSPPM;
        const ScenePtr& mScene;
        int mCurrentIteration;
        const SampleRange& mSampleRange;
        const PermutedHalton& mHalton;
        vector<uint64_t> mHaltonStartID;
        RNG mRNG;
    };

    void RayTraceTask::run(TLSPtr& tls) {
        if (mScene->getLights().empty()) {
            return;
        }
        RayTraceTLS* rayTraceTLS = static_cast<RayTraceTLS*>(tls.get());
        size_t pixelOffset = 0;
        for (int y = mSampleRange.yStart; y < mSampleRange.yEnd; ++y) {
            for (int x = mSampleRange.xStart; x  < mSampleRange.xEnd; ++x) {
                uint64_t id = mHaltonStartID[pixelOffset] + mCurrentIteration;
                mHalton.sample(&rayTraceTLS->mSample, x, y, id, &mRNG);
                mSPPM->rayTracePass(mScene, rayTraceTLS->mSample, x, y);
                pixelOffset++;
            }
        }
    }

    // Each thread owns a list of PhotonCache during photon tracing pass to
    // store the photon intersection result. The info PhotonCache carries
    // will be merged to PixelData after photon tracing pass finishes
    // (per iteration)
    struct PhotonCache {
        PhotonCache(): Phi(0.0f), Mi(0) {}
        Color Phi;
        size_t Mi;
    };


    class PhotonTraceTLS : public ThreadLocalStorage {
    public:
        PhotonTraceTLS(const SampleQuota& sampleQuota, size_t threadID,
            vector<PhotonCache>& photonCache):
            mThreadID(threadID), mPhotonCache(photonCache),
            mEmittedPhotons(0) {
            mSample.allocateQuota(sampleQuota);
        }
        Sample mSample;
        size_t mThreadID;
        vector<PhotonCache>& mPhotonCache;
        uint64_t mEmittedPhotons;
    };


    class PhotonTraceTLSManager : public TLSManager {
    public:
        PhotonTraceTLSManager(const SampleQuota& sampleQuota,
            vector<PixelData>& pixelData,
            vector<vector<PhotonCache> >& photonCache,
            uint64_t* emittedPhotons):
            mSampleQuota(sampleQuota), mNextThreadID(0), mPixelData(pixelData),
            mPhotonCache(photonCache), mEmittedPhotons(emittedPhotons) {}

        void initialize(TLSPtr& tlsPtr) {
            {
                // TODO this can probably be replaced with an atomic
                boost::lock_guard<boost::mutex> lk(mSyncTLSMutex);
                size_t threadID = mNextThreadID++;
                tlsPtr.reset(new PhotonTraceTLS(mSampleQuota, threadID,
                    mPhotonCache[threadID]));
            }
        }

        void finalize(TLSPtr& tlsPtr) {
            {
                boost::lock_guard<boost::mutex> lk(mSyncTLSMutex);
                PhotonTraceTLS* photonTraceTLS =
                    static_cast<PhotonTraceTLS*>(tlsPtr.get());
                vector<PhotonCache>& photonCache =
                    photonTraceTLS->mPhotonCache;
                for (size_t i = 0; i < mPixelData.size(); ++i) {
                    if (mPixelData[i].throughput == Color::Black) {
                        continue;
                    }
                    mPixelData[i].Phi += photonCache[i].Phi;
                    photonCache[i].Phi = Color::Black;
                    mPixelData[i].Mi += photonCache[i].Mi;
                    photonCache[i].Mi = 0;
                }
                *mEmittedPhotons += photonTraceTLS->mEmittedPhotons;
            }
        }

    private:
        boost::mutex mSyncTLSMutex;
        const SampleQuota& mSampleQuota;
        size_t mNextThreadID;
        vector<PixelData>& mPixelData;
        vector<vector<PhotonCache> >& mPhotonCache;
        uint64_t* mEmittedPhotons;
    };


    class PhotonTraceTask : public Task {
    public:
        PhotonTraceTask(SPPM* sppm, const ScenePtr& scene,
            const PermutedHalton& halton, uint64_t haltonOffset,
            uint64_t sampleNum):
            mSPPM(sppm), mScene(scene), mHalton(halton),
            mIterationOffset(0), mHaltonOffset(haltonOffset),
            mSampleNum(sampleNum) {}

        void run(TLSPtr& tls);

        void setIterationOffset(uint64_t offset) {
            mIterationOffset = offset;
        }
    private:
        uint64_t getHaltonStartID() const {
            return mIterationOffset + mHaltonOffset;
        }

    private:
        SPPM* mSPPM;
        const ScenePtr& mScene;
        const PermutedHalton& mHalton;
        uint64_t mIterationOffset;
        uint64_t mHaltonOffset;
        uint64_t mSampleNum;
        RNG mRNG;
    };

    void PhotonTraceTask::run(TLSPtr& tls) {
        if (mScene->getLights().empty()) {
            return;
        }
        PhotonTraceTLS* photonTraceTLS =
            static_cast<PhotonTraceTLS*>(tls.get());
        for (uint64_t i = 0; i < mSampleNum; ++i) {
            uint64_t id = getHaltonStartID() + i;
            mHalton.sample(&photonTraceTLS->mSample, id, &mRNG);
            mSPPM->photonTracePass(mScene, photonTraceTLS->mSample,
                photonTraceTLS->mPhotonCache);
        }
        photonTraceTLS->mEmittedPhotons += mSampleNum;
    }


    class SpatialHashGrids {
    public:
        SpatialHashGrids(const ImageRect& filmRect);

        void rebuild(vector<PixelData>& pixelData);

        bool worldToGrid(const Vector3& p,
            uint32_t* x, uint32_t* y, uint32_t* z) const;

        const vector<PixelData*>* getGrid(const Vector3& p) const;

    private:
        uint32_t hash(uint32_t x, uint32_t y, uint32_t z) const;

    private:
        vector<vector<PixelData*> > mGrids;
        BBox mGridsBBox;
        int mXRes;
        int mYRes;
        size_t mHashSize;
        float mInvGridLength;

    };

    SpatialHashGrids::SpatialHashGrids(const ImageRect& filmRect):
        mGrids(filmRect.pixelNum()),
        mXRes(filmRect.xCount), mYRes(filmRect.yCount),
        mHashSize(filmRect.pixelNum()), mInvGridLength(1.0f) {}

    void SpatialHashGrids::rebuild(vector<PixelData>& pixelData) {
        for (size_t i = 0; i < mGrids.size(); ++i) {
            vector<PixelData*>().swap(mGrids[i]);
        }
        BBox gridsBBox;
        float maxRadius = PixelData::sInvalidRadius;
        for (size_t i = 0; i < pixelData.size(); ++i) {
            if (pixelData[i].throughput != Color::Black) {
                gridsBBox.expand(pixelData[i].fragment.getPosition());
                if (pixelData[i].Ri > maxRadius) {
                    maxRadius = pixelData[i].Ri;
                }
            }
        }
        // figure out an initial radius with heuristic method
        if (isEqual(maxRadius, PixelData::sInvalidRadius)) {
            Vector3 axis = gridsBBox.pMax - gridsBBox.pMin;
            maxRadius = ((axis.x + axis.y + axis.z) / 3.0f) /
                ((mXRes + mYRes) / 2.0f) * 2.0f;
            // we did what we can, just fall back to a inital radius
            if (maxRadius == 0.0f) {
                maxRadius = 1e-5f;
            }
            // initalize PixelData with this inital radius
            for (size_t i = 0; i < pixelData.size(); ++i) {
                pixelData[i].Ri = maxRadius;
            }
        }
        // make each grid cell two times larger than max radius
        gridsBBox.expand(maxRadius);
        mGridsBBox = gridsBBox;
        mInvGridLength = 1.0f / (2.0f * maxRadius);
        // hash the input PixelData into grids
        for (size_t i = 0; i < pixelData.size(); ++i) {
            if (pixelData[i].throughput == Color::Black) {
                continue;
            }
            Vector3 r(maxRadius, maxRadius, maxRadius);
            Vector3 pMin = pixelData[i].fragment.getPosition() - r;
            Vector3 pMax = pixelData[i].fragment.getPosition() + r;
            uint32_t xMin, yMin, zMin, xMax, yMax, zMax;
            worldToGrid(pMin, &xMin, &yMin, &zMin);
            worldToGrid(pMax, &xMax, &yMax, &zMax);
            for (uint32_t z = zMin; z <= zMax; ++z) {
                for (uint32_t y = yMin; y <= yMax; ++y) {
                    for (uint32_t x = xMin; x <= xMax; ++x) {
                        uint32_t id = hash(x, y, z);
                        mGrids[id].push_back(&pixelData[i]);
                    }
                }
            }
        }
    }

    bool SpatialHashGrids::worldToGrid(const Vector3& p,
        uint32_t* x, uint32_t* y, uint32_t* z) const {
        Vector3 gridCoord = (p - mGridsBBox.pMin) * mInvGridLength;
        *x = floorInt(gridCoord.x);
        *y = floorInt(gridCoord.y);
        *z = floorInt(gridCoord.z);
        return (*x >= 0) && (*y >= 0) && (*z >= 0);
    }

    const vector<PixelData*>* SpatialHashGrids::getGrid(
        const Vector3& p) const {
        uint32_t x, y, z;
        if (!worldToGrid(p, &x, &y, &z)) {
            return NULL;
        }
        return &mGrids[hash(x, y, z)];
    }

    // spatial hash function that taken from small ppm
    // the original paper reference:
    // Optimized Spatial Hashing for Collision Detection of Deformable Objects
    // Mathias Teschner
    uint32_t SpatialHashGrids::hash(uint32_t x, uint32_t y, uint32_t z) const {
        return ((x * 73856093)^(y * 19349663)^(z * 83492791)) % mHashSize;
    }


    SPPM::SPPM(int samplePerPixel, int threadNum, int maxPathLength,
        float initialRadius):
        Renderer(samplePerPixel, threadNum),
        mMaxPathLength(maxPathLength), mHashGrids(NULL),
        mInitialRadius(initialRadius) {
        if (mInitialRadius <= 0.0f &&
            !isEqual(mInitialRadius, PixelData::sInvalidRadius)) {
            mInitialRadius = PixelData::sInvalidRadius;
        }
    }

    SPPM::~SPPM() { delete mHashGrids; mHashGrids = NULL; }

    Color SPPM::Li(const ScenePtr& scene, const RayDifferential& ray,
        const Sample& sample, const RNG& rng,
        RenderingTLS* tls) const {
        // sppm doesn't use this camera based method actually
        return Color::Black;
    }

    void SPPM::rayTracePass(const ScenePtr& scene, const Sample& sample,
        int pixelX, int pixelY) {
        // rayTracePass only needs to find pixel landing point and
        // calculate direct lighting. We'll let photon pass to handle GI
        const CameraPtr camera = scene->getCamera();
        Film* film = camera->getFilm();
        ImageRect filmRect;
        film->getImageRect(filmRect);
        size_t pOffset = filmRect.pixelToOffset(pixelX, pixelY);
        RayDifferential ray;
        camera->generateRay(sample, &ray);
        int pathLength = 0;
        Color throughput(1.0f);
        float epsilon;
        Intersection isect;
        while (pathLength < mMaxPathLength) {
            if (!scene->intersect(ray, &epsilon, &isect)) {
                // get image based lighting if ray didn't hit anything
                if (pathLength == 0) {
                    mPixelData[pOffset].Ld += throughput *
                        scene->evalEnvironmentLight(ray);
                }
                break;
            }
            pathLength++;
            // eye ray hit the light
            if (pathLength == 1) {
                mPixelData[pOffset].Ld += throughput * isect.Le(-ray.d);
            }
            isect.computeUVDifferential(ray);
            // calculate direct lighting
            LightSample ls(sample, mLightSampleIndexes[0], 0);
            BSDFSample bs(sample, mBSDFSampleIndexes[pathLength], 0);
            float pickSample =
                sample.u1D[mPickLightSampleIndexes[0].offset][0];
            mPixelData[pOffset].Ld += throughput *
                singleSampleLd(scene, ray, epsilon, isect, sample,
                ls, bs, pickSample);
            const Fragment& fragment = isect.fragment;
            const MaterialPtr& material = isect.getMaterial();
            Vector3 wo = -normalize(ray.d);
            bool isDiffuse = (material->getType() & BSDFDiffuse) != 0;
            if (isDiffuse || pathLength == mMaxPathLength - 1) {
                mPixelData[pOffset].fragment = fragment;
                mPixelData[pOffset].wo = wo;
                mPixelData[pOffset].material = material.get();
                mPixelData[pOffset].throughput = throughput;
                mPixelData[pOffset].pathLength = pathLength;
                break;
            }
            // spawn new ray based on bsdf sample and update throughput
            Vector3 wi;
            float bsdfPdf;
            BSDFType sampledType;
            Color f = material->sampleBSDF(fragment, wo, bs, &wi, &bsdfPdf,
                BSDFAll, &sampledType);
            if (f == Color::Black || bsdfPdf == 0.0f) {
                break;
            }
            throughput *= f * absdot(wi, fragment.getNormal()) / bsdfPdf;
            ray = RayDifferential(fragment.getPosition(), wi, epsilon);
        }
    }

    void SPPM::photonTracePass(const ScenePtr& scene, const Sample&sample,
        vector<PhotonCache>& photonCache) {
        // pick up a light
        float pickSample = sample.u1D[mPickLightSampleIndexes[0].offset][0];
        float pickLightPdf;
        const Light* light = scene->sampleLight(pickSample, &pickLightPdf);
        if (light == NULL ||pickLightPdf == 0.0f) {
            return;
        }
        // sample a photon from light (position, direction)
        // TODO cosolidate samplePosition/sampleDirection/eval into one
        // samplePhoton (something like that...) utility function
        LightSample ls(sample, mLightSampleIndexes[0], 0);
        Vector3 nLight;
        float pdfLightArea;
        Vector3 pLight = light->samplePosition(scene, ls, &nLight,
            &pdfLightArea);
        float pdfLightDirection;
        BSDFSample bs(sample, mBSDFSampleIndexes[0], 0);
        Vector3 dirLight = light->sampleDirection(
            nLight, bs.uDirection[0], bs.uDirection[1], &pdfLightDirection);
        // the weight of photon
        float cosTheta = light->isDelta() ? 1.0f : absdot(nLight, dirLight);
        Color photonWeight(light->eval(pLight, nLight, dirLight) * cosTheta /
            (pdfLightArea * pdfLightDirection * pickLightPdf));
        Ray ray(pLight, dirLight, 1e-5f);
        int pathLength = 0;
        float epsilon;
        Intersection isect;
        while (pathLength < mMaxPathLength) {
            if (!scene->intersect(ray, &epsilon, &isect)) {
                break;
            }
            pathLength++;
            const Fragment& fragment = isect.fragment;
            Vector3 wi = -normalize(ray.d);
            const Vector3& p = fragment.getPosition();
            if (pathLength > 1) {
                // hash p to corresponding grid
                const vector<PixelData*>* grid = mHashGrids->getGrid(p);
                if (grid != NULL) {
                    for (size_t i = 0; i < grid->size(); ++i) {
                        const PixelData* pixelData = (*grid)[i];
                        const Vector3& pixelPos =
                            pixelData->fragment.getPosition();
                        float Ri = pixelData->Ri;
                        int combinedPathLength =
                            pixelData->pathLength + pathLength;
                        if ((pixelPos - p).squaredLength() <= Ri * Ri &&
                             combinedPathLength <= mMaxPathLength) {
                             photonCache[pixelData->pixelIndex].Mi++;
                             Color fs = pixelData->material->bsdf(
                                 pixelData->fragment, pixelData->wo, wi);
                             photonCache[pixelData->pixelIndex].Phi +=
                                 fs * photonWeight;
                        }
                    }
                }
            }
            // spawn new ray based on bsdf sample and update throughput
            const MaterialPtr& material = isect.getMaterial();
            Vector3 wo;
            float bsdfPdf;
            BSDFSample bs(sample, mBSDFSampleIndexes[pathLength], 0);
            Color f = material->sampleBSDF(fragment, wi, bs, &wo, &bsdfPdf);
            if (f == Color::Black || bsdfPdf == 0.0f) {
                break;
            }
            photonWeight *= f * absdot(wo, fragment.getNormal()) / bsdfPdf;
            ray = Ray(p, wo, epsilon);
        }
    }

    void SPPM::render(const ScenePtr& scene) {
        const CameraPtr camera = scene->getCamera();
        Film* film = camera->getFilm();
        SampleQuota sampleQuota;
        querySampleQuota(scene, &sampleQuota);

        // split up film for RayTraceTask to work on
        ImageRect filmRect;
        film->getImageRect(filmRect);
        int xStart = filmRect.xStart;
        int xEnd = xStart + filmRect.xCount;
        int yStart = filmRect.yStart;
        int yEnd = yStart + filmRect.yCount;
        int tileWidth = 64;
        vector<SampleRange> sampleRanges;
        for (int y = yStart; y < yEnd; y += tileWidth) {
            for (int x = xStart; x < xEnd; x += tileWidth) {
                sampleRanges.push_back(SampleRange(
                    x, min(x + tileWidth, xEnd),
                    y, min(y + tileWidth, yEnd)));
            }
        }
        RNG rng;
        // init HaltonSampler for RayTraceTask
        PermutedHalton rayTraceHalton(sampleQuota.getDimension(), &rng);
        // init PixelData for each pixel
        mPixelData.resize(filmRect.pixelNum());
        for (size_t i = 0; i < mPixelData.size(); ++i) {
            mPixelData[i].pixelIndex = i;
            mPixelData[i].Ri = mInitialRadius;
        }
        if (mHashGrids) {
            delete mHashGrids;
        }
        mHashGrids = new SpatialHashGrids(filmRect);
        // init RayTraceTask
        vector<Task*> rayTraceTasks(sampleRanges.size());
        for (size_t i = 0; i < rayTraceTasks.size(); ++i) {
            rayTraceTasks[i] = new RayTraceTask(
                this, scene, sampleRanges[i], rayTraceHalton);
        }
        // init HaltonSampler for PhotonTraceTask
        PermutedHalton photonTraceHalton(sampleQuota.getDimension(), &rng);
        // init PhotonTraceTask
        vector<Task*> photonTraceTasks(mThreadNum);
        size_t taskPhotonSamples = max(filmRect.pixelNum() / mThreadNum, 1);
        for (size_t i = 0 ; i < photonTraceTasks.size(); ++i) {
            photonTraceTasks[i] = new PhotonTraceTask(
                this, scene, photonTraceHalton,
                i * taskPhotonSamples, taskPhotonSamples);
        }
        vector<vector<PhotonCache> > photonChaches(mThreadNum);
        for (size_t i = 0; i < photonChaches.size(); ++i) {
            photonChaches[i].resize(filmRect.pixelNum());
        }
        uint64_t emittedPhotons = 0;
        int iterationCount = mSamplePerPixel;
        for (int i = 0; i < iterationCount; ++i) {
            // ray trace pass
            RayTraceTLSManager rayTraceTLSManager(sampleQuota);
            ThreadPool rayTraceThreadPool(mThreadNum, &rayTraceTLSManager);
            rayTraceThreadPool.enqueue(rayTraceTasks);
            rayTraceThreadPool.waitForAll();
            for (size_t j = 0; j < rayTraceTasks.size(); ++j) {
                RayTraceTask* task =
                    static_cast<RayTraceTask*>(rayTraceTasks[j]);
                task->nextIteration();
            }
            // deposit visible pixels into hash grids
            mHashGrids->rebuild(mPixelData);

            // photon trace pass
            PhotonTraceTLSManager photonTraceTLSManager(sampleQuota,
                mPixelData, photonChaches, &emittedPhotons);
            ThreadPool photonTraceThreadPool(mThreadNum,
                &photonTraceTLSManager);
            photonTraceThreadPool.enqueue(photonTraceTasks);
            photonTraceThreadPool.waitForAll();
            for (size_t j = 0; j < photonTraceTasks.size(); ++j) {
                PhotonTraceTask* task =
                    static_cast<PhotonTraceTask*>(photonTraceTasks[j]);
                task->setIterationOffset(emittedPhotons);
            }
            // update sppm data (Tau, radius, Ni)
            for (size_t j = 0; j < mPixelData.size(); ++j) {
                const float alpha = 0.7f;
                if (mPixelData[j].throughput == Color::Black) {
                    continue;
                }
                if (mPixelData[j].Mi > 0) {
                    float newNi = mPixelData[j].Ni + alpha * mPixelData[j].Mi;
                    float Ri = mPixelData[j].Ri;
                    float newRi = Ri *
                        sqrtf(newNi / (mPixelData[j].Ni + mPixelData[j].Mi));
                    const Color& Tau = mPixelData[j].Tau;
                    Color newTau =
                        (Tau + mPixelData[j].throughput * mPixelData[j].Phi) *
                        (newRi / Ri) * (newRi / Ri);
                    mPixelData[j].Ni = newNi;
                    mPixelData[j].Ri = newRi;
                    mPixelData[j].Tau = newTau;
                }
                mPixelData[j].reset();
            }

            // report progress
            std::cout << "\rIteration: " << i + 1 << "/" << iterationCount;
            std::cout.flush();
            if(i == iterationCount - 1) {
                std::cout << "\rRender Complete!         " << std::endl;
                std::cout.flush();
            }
        }
        // clean up
        for (size_t i = 0; i < rayTraceTasks.size(); ++i) {
            delete rayTraceTasks[i];
        }
        rayTraceTasks.clear();
        for (size_t i = 0; i <photonTraceTasks.size(); ++i) {
            delete photonTraceTasks[i];
        }

        ImageTile tile(filmRect, film->getFilterTable());
        float invIterationCount = 1.0f / (float)iterationCount;
        for (size_t i = 0; i < mPixelData.size(); ++i) {
            int x, y;
            filmRect.offsetToPixel(i, &x, &y);
            // direct lighting from ray trace pass
            Color Ld = mPixelData[i].Ld * invIterationCount;
            float r = mPixelData[i].Ri;
            // indiret lighting from photon trace pass
            Color Lbounce = mPixelData[i].Tau / (emittedPhotons * PI * r * r);
            tile.addSample((float)x, (float)y, Ld + Lbounce);
        }

        film->mergeTile(tile);
        film->writeImage();
    }

    void SPPM::querySampleQuota(const ScenePtr& scene,
        SampleQuota* sampleQuota){

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
        // for direct lighting sampling
        mLightSampleIndexes = new LightSampleIndex[1];
        mLightSampleIndexes[0] = LightSampleIndex(sampleQuota, 1);
        // for pick up light in direct lighting
        mPickLightSampleIndexes = new SampleIndex[1];
        mPickLightSampleIndexes[0] = sampleQuota->requestOneDQuota(1);
        // for spawning new ray
        mBSDFSampleIndexes = new BSDFSampleIndex[mMaxPathLength + 1];
        for(int i = 0; i < mMaxPathLength + 1; ++i) {
            mBSDFSampleIndexes[i] = BSDFSampleIndex(sampleQuota, 1);
        }
    }

    Renderer* SPPMCreator::create(const ParamSet& params) const {
        int samplePerPixel = params.getInt("sample_per_pixel", 1);
        int threadNum = params.getInt("thread_num",
            boost::thread::hardware_concurrency());
        int maxPathLength = params.getInt("max_path_length", 5);
        float initialRadius = params.getFloat("initial_radius",
            PixelData::sInvalidRadius);
        return new SPPM(samplePerPixel, threadNum, maxPathLength,
            initialRadius);
    }
}
