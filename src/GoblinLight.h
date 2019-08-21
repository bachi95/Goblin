#ifndef GOBLIN_LIGHT_H
#define GOBLIN_LIGHT_H

#include "GoblinColor.h"
#include "GoblinGeometry.h"
#include "GoblinMaterial.h"
#include "GoblinParamSet.h"
#include "GoblinTexture.h"
#include "GoblinTransform.h"
#include "GoblinUtils.h"
#include "GoblinVector.h"

namespace Goblin {
class Ray;
class Quaternion;
class CDF1D;
class CDF2D;
class SampleQuota;
class Sample;
struct SampleIndex;

struct LightSampleIndex {
    LightSampleIndex() {}
    LightSampleIndex(SampleQuota* sampleQuota, int requestNum);
    uint32_t samplesNum;
    uint32_t componentIndex;
    uint32_t geometryIndex;
};

struct LightSample {
    LightSample(const RNG& rng);
    LightSample(const Sample& sample,
        const LightSampleIndex& index, uint32_t n);
    float uComponent;
    float uGeometry[2];
};

struct BSSRDFSampleIndex {
    BSSRDFSampleIndex() {}
    BSSRDFSampleIndex(SampleQuota* sampleQuota, int requestNum);
    LightSampleIndex lsIndex;
    uint32_t pickLightIndex;
    uint32_t pickAxisIndex;
    uint32_t discSampleIndex;
    uint32_t singleScatterIndex;
    uint32_t samplesNum;
};

struct BSSRDFSample {
    BSSRDFSample(const RNG& rng);
    BSSRDFSample(const Sample& sample,
        const BSSRDFSampleIndex& index, uint32_t n);
    float uPickLight;
    float uPickAxis;
    LightSample ls;
    float uDisc[2];
    float uSingleScatter;
};

class Light {
public:
    enum Type {
        Point = 0,
        Directional = 1,
        Spot = 2,
        Area = 3,
        IBL = 4
    };

    Light();

    virtual ~Light() {};

    // this one is only usable for IBL for now...
    // when ray doesn't intersect scene
    virtual Color Le(const Ray& ray) const { return Color::Black; }

    virtual Color sampleL(const Vector3& p, float epsilon,
        const LightSample& lightSample,
        Vector3* wi, float* pdf, Ray* shadowRay) const = 0;

    // sample a position in light surface
    virtual Vector3 samplePosition(const ScenePtr& scene,
        const LightSample& ls, Vector3* surfaceNormal,
        float* pdfArea) const = 0;

    // sample a direction with the given position in light surface
    virtual Vector3 sampleDirection(const Vector3& pSurface,
        float u1, float u2, float* pdfW) const = 0;

    // get the pdf (measured in area) this particular point
    // on light is sampled
    virtual float pdfPosition(const ScenePtr& scene,
        const Vector3& p) const = 0;

    // get the conditional pdf (measured in solid andgle)
    // this light emit particle in direction wo (under the
    // condition that it is already sampled in point p with normal n)
    virtual float pdfDirection(const Vector3& p, const Vector3& n,
        const Vector3& wo) const = 0;

    // evaluate the flux differential with regards to spatial component
    // and directional component (in area light case, this is the radiance
    // emitted on specified light surface point/direction)
    virtual Color eval(const Vector3& p, const Vector3& n,
        const Vector3& wo) const = 0;

    // get the pdf (measured in solid angle) this light get sampled
    // from point p with direction wi
    virtual float pdf(const Vector3& p, const Vector3& wi) const {
        return 0.0f;
    }

    // whether the light source is in delta distribution
    // (For example, traditional cg light like point/directional/spot)
    virtual bool isDelta() const { return true; }

    // whether the light source is in infinite distance
    // (For example, IBL and directional light)
    virtual bool isInfinite() const { return false; }

    virtual Color power(const Scene& scene) const = 0;

    virtual uint32_t getSamplesNum() const { return 1; }

    size_t getId() const { return mLightId; }

    const ParamSet& getParams() const { return mParams; }

protected:
    void setOrientation(const Vector3& dir);

protected:
    static size_t nextLightId;
    size_t mLightId;
    ParamSet mParams;
    Transform mToWorld;
};


class PointLight : public Light {
public:
    PointLight(const Color& intensity, const Vector3& position);

    Color sampleL(const Vector3& p, float epsilon,
        const LightSample& lightSample,
        Vector3* wi, float* pdf, Ray* shadowRay) const;

    Vector3 samplePosition(const ScenePtr& scene,
        const LightSample& ls, Vector3* surfaceNormal,
        float* pdfArea) const;

    Vector3 sampleDirection(const Vector3& pSurface,
        float u1, float u2, float* pdfW) const;

    float pdfPosition(const ScenePtr& scene,
        const Vector3& p) const;

    float pdfDirection(const Vector3& p, const Vector3& n,
        const Vector3& wo) const;

    Color eval(const Vector3& p, const Vector3& n,
        const Vector3& wo) const;

    Color power(const Scene& scene) const;
private:
    Color mIntensity;
};


class DirectionalLight : public Light {
public:
    DirectionalLight(const Color& radiance, const Vector3& direction);

    Color sampleL(const Vector3& p, float epsilon,
        const LightSample& lightSample,
        Vector3* wi, float* pdf, Ray* shadowRay) const;

    Vector3 samplePosition(const ScenePtr& scene,
        const LightSample& ls, Vector3* surfaceNormal,
        float* pdfArea) const;

    Vector3 sampleDirection(const Vector3& pSurface,
        float u1, float u2, float* pdfW) const;

    float pdfPosition(const ScenePtr& scene,
        const Vector3& p) const;

    float pdfDirection(const Vector3& p, const Vector3& n,
        const Vector3& wo) const;

    Color eval(const Vector3& p, const Vector3& n,
        const Vector3& wo) const;

    bool isInfinite() const { return true; }

    Color power(const Scene& scene) const;

    Vector3 getDirection() const { return mToWorld.onVector(Vector3::UnitZ); }

    void setDirection(const Vector3& dir) { setOrientation(dir); }

private:
    Color mRadiance;
};


class SpotLight : public Light {
public:
    SpotLight(const Color& intensity, const Vector3& position,
        const Vector3& dir, float cosThetaMax, float cosFalloffStart);

    Color sampleL(const Vector3& p, float epsilon,
        const LightSample& lightSample,
        Vector3* wi, float* pdf, Ray* shadowRay) const;

    Vector3 samplePosition(const ScenePtr& scene,
        const LightSample& ls, Vector3* surfaceNormal,
        float* pdfArea) const;

    Vector3 sampleDirection(const Vector3& pSurface,
        float u1, float u2, float* pdfW) const;

    float pdfPosition(const ScenePtr& scene,
        const Vector3& p) const;

    float pdfDirection(const Vector3& p, const Vector3& n,
        const Vector3& wo) const;

    Color eval(const Vector3& p, const Vector3& n,
        const Vector3& wo) const;

    Color power(const Scene& scene) const;
private:
    float falloff(const Vector3& w) const;
private:
    Color mIntensity;
    float mCosThetaMax;
    float mCosFalloffStart;
};


class GeometrySet {
public:
    GeometrySet(const Geometry* geometry);

    ~GeometrySet();

    Vector3 sample(const Vector3& p, const LightSample& lightSample,
        Vector3* normal) const;

    Vector3 sample(const LightSample& lightSample,
        Vector3* normal) const;

    float pdf(const Vector3& p, const Vector3& wi) const;

    float area() const { return mSumArea; }

private:
    GeometryList mGeometries;
    std::vector<float> mGeometriesArea;
    float mSumArea;
    CDF1D* mAreaDistribution;
};


class AreaLight : public Light {
public:
    AreaLight(const Color& Le, const Geometry* geometry,
        const Transform& toWorld, uint32_t samplesNum);

    ~AreaLight();

    Color sampleL(const Vector3& p, float epsilon,
        const LightSample& lightSample,
        Vector3* wi, float* pdf, Ray* shadowRay) const;

    Vector3 samplePosition(const ScenePtr& scene,
        const LightSample& ls, Vector3* surfaceNormal,
        float* pdfArea) const;

    Vector3 sampleDirection(const Vector3& surfaceNormal,
        float u1, float u2, float* pdfW) const;

    float pdfPosition(const ScenePtr& scene,
        const Vector3& p) const;

    float pdfDirection(const Vector3& p, const Vector3& n,
        const Vector3& wo) const;

    Color eval(const Vector3& p, const Vector3& n,
        const Vector3& wo) const;

    float pdf(const Vector3& p, const Vector3& wi) const;

    bool isDelta() const { return false; }
    /*
        * ps: point on area light surface
        * ns: normal on area light surface
        * w: direction L leaving surface
        */
    Color L(const Vector3& ps, const Vector3& ns, const Vector3& w) const;

    Color power(const Scene& scene) const;

    uint32_t getSamplesNum() const { return mSamplesNum; }
private:
    Color mLe;
    uint32_t mSamplesNum;
    GeometrySet* mGeometrySet;
};


class ImageBasedLight : public Light {
public:
    ImageBasedLight(const std::string& radianceMap, const Color& filter,
        const Quaternion& orientation = Quaternion::Identity,
        uint32_t samplesNum = 1);

    ~ImageBasedLight();

    Color Le(const Ray& ray) const;

    Color sampleL(const Vector3& p, float epsilon,
        const LightSample& lightSample,
        Vector3* wi, float* pdf, Ray* shadowRay) const;

    Vector3 samplePosition(const ScenePtr& scene,
        const LightSample& ls, Vector3* surfaceNormal,
        float* pdfArea) const;

    Vector3 sampleDirection(const Vector3& pSurface,
        float u1, float u2, float* pdfW) const;

    float pdfPosition(const ScenePtr& scene,
        const Vector3& p) const;

    float pdfDirection(const Vector3& p, const Vector3& n,
        const Vector3& wo) const;

    Color eval(const Vector3& p, const Vector3& n,
        const Vector3& wo) const;

    float pdf(const Vector3& p, const Vector3& wi) const;

    bool isDelta() const { return false; }

    bool isInfinite() const { return true; }

    Color power(const Scene& scene) const;

    uint32_t getSamplesNum() const { return mSamplesNum; }
private:
    MIPMap<Color>* mRadiance;
    CDF2D* mDistribution;
    Color mAverageRadiance;
    uint32_t mSamplesNum;
    int mSampleMIPLevel;
};


Light* createPointLight(const ParamSet& params,
    const SceneCache& sceneCache);

Light* createDirectionalLight(const ParamSet& params,
    const SceneCache& sceneCache);

Light* createSpotLight(const ParamSet& params,
	const SceneCache& sceneCache);

Light* createAreaLight(const ParamSet& params,
    const SceneCache& sceneCache);

Light* createImageBasedLight(const ParamSet& params,
    const SceneCache& sceneCache);

}

#endif //GOBLIN_LIGHT_H