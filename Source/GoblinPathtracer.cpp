#include "GoblinPathtracer.h"
#include "GoblinRay.h"

namespace Goblin {
    static bool isOpaque(const Primitive* p, const Ray& r) {
        return !(p->getMaterial()->getType() & BSDFNull);
    }
    
    static bool notOpaque(const Primitive* p, const Ray& r) {
        return !isOpaque(p, r);
    }

    PathTracer::PathTracer(int samplePerPixel, int threadNum, 
        int maxRayDepth, int bssrdfSampleNum): 
        Renderer(samplePerPixel, threadNum),
        mMaxRayDepth(maxRayDepth),
        mBssrdfSampleNum(bssrdfSampleNum) {}

    PathTracer::~PathTracer() {}

    Color PathTracer::evalAttenuation(const ScenePtr& scene, 
        const Ray& ray, const BSDFSample& bs) const {
        Color throughput(1.0f);
        float maxt = ray.maxt;
        Ray currentRay(ray);
        float epsilon;
        Intersection intersection;
        while(true) {
            if (!scene->intersect(currentRay,
                &epsilon, &intersection, &notOpaque)) {
                break;
            } 
            const MaterialPtr& material = intersection.getMaterial();
            const Fragment& fragment = intersection.fragment;
            Vector3 wo = -currentRay.d;   
            Vector3 wi;
            float pdf;
            throughput *= material->sampleBSDF(fragment, wo, bs, 
                &wi, &pdf, BSDFNull);
            if (throughput == Color::Black) {
                break;
            }
            // move forward ray's start point
            currentRay.mint = currentRay.maxt + epsilon;
            currentRay.maxt = maxt;
        }
        return throughput;
    }

    Color PathTracer::Li(const ScenePtr& scene, const RayDifferential& ray, 
        const Sample& sample, const RNG& rng,
        RenderingTLS* tls) const {
        const vector<Light*>& lights = scene->getLights();
        if (lights.size() == 0) {
            return Color(0.0f);
        }

        Color Li(0.0f);
        float epsilon;
        Intersection intersection;
        if (!scene->intersect(ray, &epsilon, &intersection)) {
            // get image based lighting if ray didn't hit anything
            Li += scene->evalEnvironmentLight(ray);
            return Li;
        }
        // if intersect an area light
        Li += intersection.Le(-ray.d);
        // subsurface scattering 
        Li += Lsubsurface(scene, intersection, -ray.d, sample, 
            &mBSSRDFSampleIndex, tls);

        RayDifferential currentRay = ray;
        vector<Ray> debugRays;
        Color throughput(1.0f);
        bool firstBounce = true;
        for (int bounces = 0; bounces < mMaxRayDepth; ++bounces) {
            intersection.computeUVDifferential(currentRay);
            LightSample ls(sample, mLightSampleIndexes[bounces], 0);
            BSDFSample bs(sample, mBSDFSampleIndexes[bounces], 0);
            float pickSample = 
                sample.u1D[mPickLightSampleIndexes[bounces].offset][0];
            float pickLightPdf;
            const Light* light = scene->sampleLight(pickSample, &pickLightPdf);
            // direct lighting
            Color Ld(0.0f);
            const MaterialPtr& material = 
                intersection.primitive->getMaterial();
            const Fragment& fragment = intersection.fragment;
            Vector3 wo = -currentRay.d;
            Vector3 wi;
            Vector3 p = fragment.getPosition();
            Vector3 n = fragment.getNormal();
            float lightPdf, bsdfPdf;
            Ray shadowRay;
            // lighting sample
            Color L = light->sampleL(p, epsilon, ls, &wi, &lightPdf, &shadowRay);
            if (L != Color::Black && lightPdf > 0.0f) {
                Color f = material->bsdf(fragment, wo, wi);
                if (f != Color::Black && 
                    !scene->intersect(shadowRay, &isOpaque)) {
                    // calculate the transmittance alone index-matched material
                    Color tr = evalAttenuation(scene, shadowRay,
                        BSDFSample(rng)); 
                    // we don't do MIS for delta distribution light
                    // since there is only one sample need for it
                    if (light->isDelta()) {
                        Ld += f * tr * L * absdot(n, wi) / lightPdf;
                    } else {
                        bsdfPdf = material->pdf(fragment, wo, wi);
                        float lWeight = powerHeuristic(1, lightPdf, 1, bsdfPdf);
                        Ld += f * tr * L * absdot(n, wi) * lWeight / lightPdf;
                    }
                }
            }
            // bsdf sample
            BSDFType sampledType;
            Color f = material->sampleBSDF(fragment, wo, bs, 
                &wi, &bsdfPdf, BSDFAll, &sampledType);
            if (f != Color::Black && bsdfPdf > 0.0f) {
                // this sample steps on an index-matched BSDF
                // should punch through it with attenuation accounted
                if (sampledType == BSDFNull) {
                    throughput *= (f / bsdfPdf);
                    currentRay = RayDifferential(p, wi, epsilon);
                    if (!scene->intersect(currentRay, &epsilon, &intersection)) {
                        // primary ray need to evaluate image based lighting
                        // in this case
                        if (firstBounce) {
                            Li += throughput * 
                                scene->evalEnvironmentLight(currentRay);
                        }
                        break;
                    }
                    continue;
                }

                // calculate the misWeight if it's not a specular material
                // otherwise we should got 0 Ld from light sample earlier,
                // and count on this part for all the Ld contribution
                float fWeight = 1.0f;
                if (!(sampledType & BSDFSpecular)) {
                    lightPdf = light->pdf(p, wi);
                    fWeight = powerHeuristic(1, bsdfPdf, 1, lightPdf);
                }
                Intersection lightIntersect;
                float lightEpsilon;
                Ray r(p, wi, epsilon);
                if (scene->intersect(r, &lightEpsilon, 
                    &lightIntersect, &isOpaque)) {
                    Color tr = evalAttenuation(scene, r, BSDFSample(rng));
                    if (lightIntersect.primitive->getAreaLight() == light) {
                        Color Li = lightIntersect.Le(-wi);
                        if (Li != Color::Black) {
                            Ld += f * tr * Li * absdot(wi, n) * fWeight / bsdfPdf;
                        }
                    }
                } else {
                    // the radiance contribution from IBL
                    Color tr = evalAttenuation(scene, r, BSDFSample(rng));
                    Ld += f * tr * light->Le(r) * fWeight / bsdfPdf;
                }
            }
            Li += throughput * Ld / pickLightPdf;

            // indirect lighting
            if ( f == Color::Black || bsdfPdf == 0.0f) {
                break;
            }
            throughput *= f * absdot(wi, n) / bsdfPdf;
            currentRay = RayDifferential(p, wi, epsilon);
            if (!scene->intersect(currentRay, &epsilon, &intersection)) {
                break;
            }
            firstBounce = false;
            //debugRays.push_back(currentRay);
        }

        return Li;        
    }

    void PathTracer::querySampleQuota(const ScenePtr& scene, 
            SampleQuota* sampleQuota) {
        if (mLightSampleIndexes) {
            delete [] mLightSampleIndexes;
            mLightSampleIndexes = NULL;
        }
        if (mBSDFSampleIndexes) {
            delete [] mBSDFSampleIndexes;
            mBSDFSampleIndexes = NULL;
        }
        if (mPickLightSampleIndexes) {
            delete [] mPickLightSampleIndexes;
            mPickLightSampleIndexes = NULL;
        }

        int bounces = mMaxRayDepth;
        mLightSampleIndexes = new LightSampleIndex[bounces];
        mBSDFSampleIndexes = new BSDFSampleIndex[bounces];
        mPickLightSampleIndexes = new SampleIndex[bounces];
        for (int i = 0; i < bounces; ++i) {
            mLightSampleIndexes[i] = LightSampleIndex(sampleQuota, 1);
            mBSDFSampleIndexes[i] = BSDFSampleIndex(sampleQuota, 1);
            mPickLightSampleIndexes[i] = sampleQuota->requestOneDQuota(1);
        }

        mBSSRDFSampleIndex = BSSRDFSampleIndex(sampleQuota, 
            mBssrdfSampleNum);
    }

    Renderer* PathTracerCreator::create(const ParamSet& params) const {
        int samplePerPixel = params.getInt("sample_per_pixel", 1);
        int threadNum = params.getInt("thread_num", getMaxThreadNum());
        int maxRayDepth = params.getInt("max_ray_depth", 5);
        int bssrdfSampleNum = params.getInt("bssrdf_sample_num", 4);
        return new PathTracer(samplePerPixel, threadNum, 
            maxRayDepth, bssrdfSampleNum);
    }

}
