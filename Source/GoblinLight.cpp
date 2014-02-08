#include "GoblinLight.h"
#include "GoblinRay.h"
#include "GoblinSampler.h"

namespace Goblin {

    LightSampleIndex::LightSampleIndex(Sampler* sampler, int requestNum) {
        SampleIndex oneDIndex = sampler->requestOneDQuota(requestNum);
        SampleIndex twoDIndex = sampler->requestTwoDQuota(requestNum);
        // theoretically this two should be the same...
        // this is just a paranoid double check
        samplesNum = min(oneDIndex.sampleNum, twoDIndex.sampleNum);
        componentIndex = oneDIndex.offset;
        geometryIndex = twoDIndex.offset;
    }

    LightSample::LightSample() {
        uComponent = randomFloat();
        uGeometry[0] = randomFloat();
        uGeometry[1] = randomFloat();
    }

    LightSample::LightSample(const Sample& sample, 
        const LightSampleIndex& index, uint32_t n) {
        uComponent = sample.u1D[index.componentIndex][n];
        uGeometry[0] = sample.u2D[index.geometryIndex][2 * n];
        uGeometry[1] = sample.u2D[index.geometryIndex][2 * n + 1];
    }

    PointLight::PointLight(const Color& I, const Vector3& P):
    intensity(I), position(P) {
        mParams.setInt("type", Point);
        mParams.setColor("intensity", intensity);
        mParams.setVector3("position", position);
    }

    Color PointLight::sampleL(const Vector3& p, float epsilon, 
        const LightSample& lightSample,
        Vector3* wi, float* pdf, Ray* shadowRay) const {
        Vector3 dir = position - p;
        *wi = normalize(dir);
        *pdf = 1.0f;
        shadowRay->o = p;
        shadowRay->d = *wi;
        shadowRay->mint = epsilon;
        float squaredDistance = squaredLength(dir);
        shadowRay->maxt = sqrt(squaredDistance) - epsilon;
        return intensity / squaredDistance;
    }

    DirectionalLight::DirectionalLight(const Color& R, const Vector3& D):
    radiance(R), direction(D) {
        mParams.setInt("type", Directional);
        mParams.setColor("radiance", radiance);
        mParams.setVector3("direction", direction);
    }

    Color DirectionalLight::sampleL(const Vector3& p, float epsilon, 
        const LightSample& lightSample,
        Vector3* wi, float* pdf, Ray* shadowRay) const {
        *wi = -direction;
        *pdf = 1.0f;
        shadowRay->o = p;
        shadowRay->d = *wi;
        shadowRay->mint = epsilon;
        return radiance;
    }

    GeometrySet::GeometrySet(const GeometryPtr& geometry):
        mAreaDistribution(NULL), mSumArea(0.0f) {
        if(geometry->intersectable()) {
            mGeometries.push_back(geometry);
        } else {
            geometry->refine(mGeometries);
        }
        mSumArea = 0.0f;
        mGeometriesArea.resize(mGeometries.size());
        for(size_t i = 0; i < mGeometries.size(); ++i) {
            float area = mGeometries[i]->area();
            mGeometriesArea[i] = area;
            mSumArea += area;
        }
        mAreaDistribution = new CDF1D(mGeometriesArea);
    }

    GeometrySet::~GeometrySet() {
        if(mAreaDistribution != NULL ) {
            delete mAreaDistribution;
            mAreaDistribution = NULL;
        }
    }

    Vector3 GeometrySet::sample(const Vector3& p, 
        const LightSample& lightSample,
        Vector3* normal) const {
        // pick up a geometry to sample based on area distribution
        float uComp = lightSample.uComponent;
        int geoIndex = mAreaDistribution->sampleDiscrete(uComp);
        // sample out ps from picked up geometry surface
        float u1 = lightSample.uGeometry[0];
        float u2 = lightSample.uGeometry[1];
        Vector3 ps = mGeometries[geoIndex]->sample(p, u1, u2, normal);
        return ps;
    }

    float GeometrySet::pdf(const Vector3& p, const Vector3& wi) const {
        float pdf = 0.0f;
        for(size_t i = 0; i < mGeometries.size(); ++i) {
            pdf += mGeometriesArea[i] * mGeometries[i]->pdf(p, wi);    
        }
        pdf /= mSumArea;
        return pdf;
    }


    AreaLight::AreaLight(const Color& Le, const GeometryPtr& geometry,
        const Transform& toWorld, uint32_t samplesNum): mLe(Le), 
        mToWorld(toWorld), mSamplesNum(samplesNum) {
        mGeometrySet = new GeometrySet(geometry);
    }

    AreaLight::~AreaLight() {
        if(mGeometrySet != NULL) {
            delete mGeometrySet;
            mGeometrySet = NULL;
        }
    }

    Color AreaLight::L(const Vector3& ps, const Vector3& ns, 
        const Vector3& w) const {
        return dot(ns, w) > 0.0f ? mLe : Color::Black;
    }

    Color AreaLight::sampleL(const Vector3& p, float epsilon,
        const LightSample& lightSample,
        Vector3* wi, float* pdf, Ray* shadowRay) const {
        // transform world space p to local space since all GeometrySet methods
        // are in local space
        Vector3 pLocal = mToWorld.invertPoint(p);
        Vector3 nsLocal;
        Vector3 psLocal = mGeometrySet->sample(pLocal, lightSample, &nsLocal);
        Vector3 wiLocal = normalize(psLocal - pLocal);
        *pdf = mGeometrySet->pdf(pLocal, wiLocal);
        // transform
        Vector3 ps = mToWorld.onPoint(psLocal); 
        Vector3 ns = normalize(mToWorld.onNormal(nsLocal));
        *wi = normalize(ps - p);

        shadowRay->o = p;
        shadowRay->d = *wi;
        shadowRay->mint = epsilon;
        shadowRay->maxt = length(ps - p) - epsilon;

        return L(ps, ns, -*wi);
    }

    float AreaLight::pdf(const Vector3& p, const Vector3& wi) const {
        Vector3 pLocal = mToWorld.invertPoint(p);
        Vector3 wiLocal = mToWorld.invertVector(wi);
        return mGeometrySet->pdf(pLocal, wiLocal);
    }
}

