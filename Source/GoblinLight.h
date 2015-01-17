#ifndef GOBLIN_LIGHT_H
#define GOBLIN_LIGHT_H

#include "GoblinColor.h"
#include "GoblinFactory.h"
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
        virtual ~Light() {};
        // this one is only usable for IBL for now...
        // when ray doesn't intersect scene
        virtual Color Le(const Ray& ray, float pdf = INFINITY, 
            BSDFType type = BSDFAll) const;
        virtual Color sampleL(const Vector3& p, float epsilon, 
            const LightSample& lightSample, 
            Vector3* wi, float* pdf, Ray* shadowRay) const = 0;
        // sample a light ray leaving the light surface
        virtual Color sampleL(const ScenePtr& scene, const LightSample& ls,
            float u1, float u2, Ray* ray, float* pdf = NULL) const = 0;
        virtual float pdf(const Vector3& p, const Vector3& wi) const;
        virtual bool isDelta() const;
        virtual Color power(const ScenePtr& scene) const = 0;
        virtual uint32_t getSamplesNum() const;
        const ParamSet& getParams() const;
    protected:
        void setOrientation(const Vector3& dir);
    protected:
        ParamSet mParams;
        Transform mToWorld;
    };

    inline Color Light::Le(const Ray& ray, float pdf, BSDFType type) const {
        return Color::Black;
    }

    inline bool Light::isDelta() const {
        return true;
    }

    inline float Light::pdf(const Vector3& p, const Vector3& wi) const {
        return 0.0f;
    }

    inline uint32_t Light::getSamplesNum() const {
        return 1;
    }

    inline const ParamSet& Light::getParams() const {
        return mParams;
    }

    class PointLight : public Light {
    public:
        PointLight(const Color& intensity, const Vector3& position);
        Color sampleL(const Vector3& p, float epsilon, 
            const LightSample& lightSample, 
            Vector3* wi, float* pdf, Ray* shadowRay) const;
        Color sampleL(const ScenePtr& scene, const LightSample& ls,
            float u1, float u2, Ray* ray, float* pdf = NULL) const;
        Color power(const ScenePtr& scene) const;
    private:
        Color mIntensity;
    };


    class DirectionalLight : public Light {
    public:
        DirectionalLight(const Color& radiance, const Vector3& direction);
        Color sampleL(const Vector3& p, float epsilon,
            const LightSample& lightSample, 
            Vector3* wi, float* pdf, Ray* shadowRay) const;
        Color sampleL(const ScenePtr& scene, const LightSample& ls,
            float u1, float u2, Ray* ray, float* pdf = NULL) const;
        Color power(const ScenePtr& scene) const;
        Vector3 getDirection() const;
        void setDirection(const Vector3& dir);
    private:
        Color mRadiance;
    };

    inline Vector3 DirectionalLight::getDirection() const {
        return mToWorld.onVector(Vector3::UnitZ);
    }

    inline void DirectionalLight::setDirection(const Vector3& dir) {
        setOrientation(dir);
    }


    class SpotLight : public Light {
    public:
        SpotLight(const Color& intensity, const Vector3& position, 
            const Vector3& dir, float cosThetaMax, float cosFalloffStart);
        Color sampleL(const Vector3& p, float epsilon, 
            const LightSample& lightSample, 
            Vector3* wi, float* pdf, Ray* shadowRay) const;
        Color sampleL(const ScenePtr& scene, const LightSample& ls,
            float u1, float u2, Ray* ray, float* pdf = NULL) const;
        Color power(const ScenePtr& scene) const;
    private:
        float falloff(const Vector3& w) const;
    private:
        Color mIntensity;
        float mCosThetaMax;
        float mCosFalloffStart;
    };


    class GeometrySet {
    public:
        GeometrySet(const GeometryPtr& geometry);
        ~GeometrySet();
        Vector3 sample(const Vector3& p, const LightSample& lightSample,
            Vector3* normal) const;
        Vector3 sample(const LightSample& lightSample, 
            Vector3* normal) const;
        float pdf(const Vector3& p, const Vector3& wi) const;
        float area() const;
    private:
        GeometryList mGeometries;
        vector<float> mGeometriesArea;
        float mSumArea;
        CDF1D* mAreaDistribution;
    };

    inline float GeometrySet::area() const {
        return mSumArea;
    }


    class AreaLight : public Light {
    public:
        AreaLight(const Color& Le, const GeometryPtr& geometry,
            const Transform& toWorld, uint32_t samplesNum);
        ~AreaLight();
        Color sampleL(const Vector3& p, float epsilon, 
            const LightSample& lightSample, 
            Vector3* wi, float* pdf, Ray* shadowRay) const;
        Color sampleL(const ScenePtr& scene, const LightSample& ls,
            float u1, float u2, Ray* ray, float* pdf = NULL) const;
        float pdf(const Vector3& p, const Vector3& wi) const;
        bool isDelta() const;
        /*
         * ps: point on area light surface
         * ns: normal on area light surface
         * w: direction L leaving surface
         */
        Color L(const Vector3& ps, const Vector3& ns, const Vector3& w) const;
        Color power(const ScenePtr& scene) const;
        uint32_t getSamplesNum() const;
    private:
        Color mLe;
        uint32_t mSamplesNum;
        GeometrySet* mGeometrySet;
    };

    inline bool AreaLight::isDelta() const {
        return false;
    }

    inline uint32_t AreaLight::getSamplesNum() const {
        return mSamplesNum;
    }


    class ImageBasedLight : public Light {
    public:
        ImageBasedLight(const string& radianceMap, const Color& filter,
            const Quaternion& orientation = Quaternion::Identity,
            int samplePerPixel = 1);
        ~ImageBasedLight();
        Color Le(const Ray& ray, float pdf = INFINITY, 
            BSDFType type = BSDFAll) const;
        Color sampleL(const Vector3& p, float epsilon,
            const LightSample& lightSample,
            Vector3* wi, float* pdf, Ray* shadowRay) const;
        Color sampleL(const ScenePtr& scene, const LightSample& ls,
            float u1, float u2, Ray* ray, float* pdf = NULL) const;
        float pdf(const Vector3& p, const Vector3& wi) const;
        bool isDelta() const;
        Color power(const ScenePtr& scene) const;
    private:
        int getLodLevel(float pdfUV) const;
    private:
        vector<ImageBuffer<Color>* > mRadiance;
        CDF2D* mDistribution;
        Color mAverageRadiance;
        int mSamplePerPixel;
    };

    inline bool ImageBasedLight::isDelta() const {
        return false;
    }


    class PointLightCreator : public 
        Creator<Light , const ParamSet&, const SceneCache&> {
    public:
        Light* create(const ParamSet& params, 
            const SceneCache& sceneCache) const;
    };


    class DirectionalLightCreator : public 
        Creator<Light , const ParamSet&, const SceneCache&> {
    public:
        Light* create(const ParamSet& params, 
            const SceneCache& sceneCache) const;
    };


    class SpotLightCreator : public 
        Creator<Light , const ParamSet&, const SceneCache&> {
    public:
        Light* create(const ParamSet& params, 
            const SceneCache& sceneCache) const;
    };


    class ImageBasedLightCreator : public 
        Creator<Light , const ParamSet&, const SceneCache&> {
    public:
        Light* create(const ParamSet& params, 
            const SceneCache& sceneCache) const;
    };

}
#endif //GOBLIN_LIGHT_H
