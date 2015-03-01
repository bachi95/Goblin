#ifndef GOBLIN_MATERIAL_H
#define GOBLIN_MATERIAL_H

#include "GoblinColor.h"
#include "GoblinFactory.h"
#include "GoblinTexture.h"
#include "GoblinUtils.h"

namespace Goblin {
    class Fragment;
    class Vector3;
    class Matrix3;
    class Ray;
    class SampleQuota;
    class Sample;
    struct BSSRDFSample;
    
    enum BSDFType {
        BSDFReflection = 1 << 0,
        BSDFTransmission = 1 << 1,
        BSDFDiffuse = 1 << 2,
        BSDFGlossy = 1 << 3,
        BSDFSpecular = 1 << 4,
        // alpha masking material
        BSDFNull = 1 << 5,
        BSDFAll = BSDFReflection | BSDFTransmission |
            BSDFDiffuse | BSDFGlossy | BSDFSpecular | BSDFNull
    };

    enum FresnelType {
        Dieletric,
        Conductor
    };

    enum BSSRDFSampleAxis {
        UAxis,
        VAxis,
        NAxis
    };

    struct BSDFSampleIndex {
        BSDFSampleIndex() {}
        BSDFSampleIndex(SampleQuota* sampleQuota, int requestNum);
        uint32_t samplesNum;
        uint32_t componentIndex;
        uint32_t directionIndex;
    };

    struct BSDFSample {
        BSDFSample(const RNG& rng);
        BSDFSample(const Sample& sample,
            const BSDFSampleIndex& index, uint32_t n);
        float uComponent;
        float uDirection[2];
    };

    class BSSRDF {
    public:
        BSSRDF(const ColorTexturePtr& absorb, 
            const ColorTexturePtr& scatterPrime, float eta, float g = 0.0f);
        BSSRDF(const Color& Kd, const Color& diffuseMeanFreePath, 
            float eta, float g = 0.0f);
        // diffusion dipole approximation part
        Color Rd(const Fragment& fragment, float d2) const;
        float MISWeight(const Fragment& fo, const Fragment& fi,
            BSSRDFSampleAxis mainAxis, float pdf, 
            float sigmaTr, float Rmax) const;
        BSSRDFSampleAxis sampleProbeRay(const Fragment& fragment,
            const BSSRDFSample& sample, float sigmaTr, float Rmax,
            Ray* probeRay, float* pdf) const;
        Color getAttenuation(const Fragment& fragment) const;
        Color getScatter(const Fragment& fragment) const;
        Color getSigmaTr(const Fragment& fragment) const;
        float getEta() const;
        float phase(const Vector3& wi, const Vector3& wo) const;
        static float Fdr(float eta);
    private:
        static void convertFromDiffuse(const Color& Kd, 
            const Color& diffuseMeanFreePath, float A, 
            Color* absorb, Color* scatterPrime);
        static float diffuseReflectance(float alphaPrime, float A);
    private:
        ColorTexturePtr mAbsorb;
        ColorTexturePtr mScatterPrime;
        float mEta;
        float mA;
        float mG;
    };

    inline float BSSRDF::Fdr(float eta) {
        // see Donner. C 2006 Chapter 5
        // the internal Fresnel reflectivity 
        // approximated with a simple polynomial expansion
        if(eta < 1.0f) {
            return -0.4399f + 0.7099f / eta - 0.3319f / (eta * eta) +
                0.0636f / (eta * eta * eta);
        } else {
            return -1.4399f / (eta * eta) + 0.7099f / eta + 0.6681f +
                0.0636f * eta;
        }
    }

    inline Color BSSRDF::getAttenuation(const Fragment& fragment) const {
        return getScatter(fragment) + mAbsorb->lookup(fragment);
    }

    inline Color BSSRDF::getScatter(const Fragment& fragment) const {
        // scatterPrime = scatter * (1 - g)
        return mScatterPrime->lookup(fragment) / (1.0f - mG);
    }

    inline float BSSRDF::getEta() const {
        return mEta;
    }

    struct BumpShaders {
        BumpShaders(const FloatTexturePtr& bump = FloatTexturePtr(), 
            const ColorTexturePtr& normal = ColorTexturePtr()):
            bumpMap(bump), normalMap(normal) {}
        void evaluate(Fragment* fragment) const;
        FloatTexturePtr bumpMap;
        ColorTexturePtr normalMap;
    };

    // This serves the purpose of what GPU shader usually do
    // perturb the geometry info(position, normal, tangent....)
    // calculate the bsdf value, hold the texture reference...etc
    class Material {
    public:
        Material(BSDFType type, const BumpShaders& bumpShaders);
        virtual ~Material() {}

        virtual void perturb(Fragment* fragment) const;

        virtual Color bsdf(const Fragment& fragment, const Vector3& wo, 
            const Vector3& wi, BSDFType type = BSDFAll) const = 0;

        virtual Color sampleBSDF(const Fragment& fragment, 
            const Vector3& wo, const BSDFSample& bsdfSample, Vector3* wi, 
            float* pdf, BSDFType type = BSDFAll, 
            BSDFType* sampledType = NULL) const = 0;

        virtual float pdf(const Fragment& fragment,
            const Vector3& wo, const Vector3& wi,
            BSDFType type = BSDFAll) const = 0;

        virtual const BSSRDF* getBSSRDF() const;

        BSDFType getType() const;

        // material util to get the fresnel factor
        static float fresnelDieletric(float cosi, float etai, float etat);

        static float fresnelConductor(float cosi, float eta, float k);

    protected:
        bool matchType(BSDFType type, BSDFType toMatch) const;

        BSDFType getSampleType(const Vector3& wo, const Vector3& wi, 
            const Fragment& f, BSDFType type) const;

        bool sameHemisphere(const Fragment& fragment,
            const Vector3& wo, const Vector3& wi) const;

        float specularReflectDieletric(const Fragment& fragment, 
            const Vector3& wo, Vector3* wi, float etai, float etat) const;

        float specularReflectConductor(const Fragment& fragment,
            const Vector3& wo, Vector3* wi, float eta, float k) const;

        float specularRefract(const Fragment& fragment, const Vector3& wo,
            Vector3* wi, float etai, float etat) const;

    protected:
        BSDFType mType;

    private:
        BumpShaders mBumpShaders;
    };

    inline Material::Material(BSDFType type, const BumpShaders& bumpShaders):
        mType(type), mBumpShaders(bumpShaders) {}

    inline bool Material::matchType(BSDFType type, BSDFType toMatch) const {
        return (type & toMatch) == toMatch;
    }

    inline const BSSRDF* Material::getBSSRDF() const {
        return NULL;
    }

    inline BSDFType Material::getType() const {
        return mType;
    }

    Vector3 specularRefract(const Vector3& wo, const Vector3& n, 
        float etai, float etat);


    typedef boost::shared_ptr<Material> MaterialPtr;

    class LambertMaterial : public Material {
    public:
        LambertMaterial(const ColorTexturePtr& Kd, 
            const BumpShaders& bumpShaders = BumpShaders());
        Color bsdf(const Fragment& fragment, const Vector3& wo, 
            const Vector3& wi, BSDFType type) const;

        Color sampleBSDF(const Fragment& fragment, 
            const Vector3& wo, const BSDFSample& bsdfSample, Vector3* wi, 
            float* pdf, BSDFType type, BSDFType* sampledType) const;

        float pdf(const Fragment& fragment, 
            const Vector3& wo, const Vector3& wi, BSDFType type) const;

    private:
        ColorTexturePtr mDiffuseFactor;
    };

    inline LambertMaterial::LambertMaterial(const ColorTexturePtr& Kd,
        const BumpShaders& bumpShaders):
        Material(BSDFType(BSDFDiffuse | BSDFReflection), bumpShaders), 
        mDiffuseFactor(Kd) {}

    class BlinnMaterial : public Material {
    public:
        BlinnMaterial(const ColorTexturePtr& Kg, const FloatTexturePtr& exp, 
            float index, const BumpShaders& bumpShaders = BumpShaders());
        BlinnMaterial(const ColorTexturePtr& Kg, const FloatTexturePtr& exp, 
            float index, float absorption, 
            const BumpShaders& bumpShaders = BumpShaders());

        Color bsdf(const Fragment& fragment, const Vector3& wo, 
            const Vector3& wi, BSDFType type) const;

        Color sampleBSDF(const Fragment& fragment, 
            const Vector3& wo, const BSDFSample& bsdfSample, Vector3* wi, 
            float* pdf, BSDFType type, BSDFType* sampledType) const;

        float pdf(const Fragment& fragment, 
            const Vector3& wo, const Vector3& wi, BSDFType type) const;
           
    private:
        ColorTexturePtr mGlossyFactor;
        FloatTexturePtr mExp;
        float mEta, mK;
        FresnelType mFresnelType;
    };

    inline BlinnMaterial::BlinnMaterial(const ColorTexturePtr& Kg, 
        const FloatTexturePtr& exponent, float index, 
        const BumpShaders& bumpShaders):
        Material(BSDFType(BSDFGlossy | BSDFReflection), bumpShaders), 
        mGlossyFactor(Kg), mExp(exponent), mEta(index),
        mFresnelType(Dieletric) {}

    inline BlinnMaterial::BlinnMaterial(const ColorTexturePtr& Kg, 
        const FloatTexturePtr& exponent, float index, float absorption,
        const BumpShaders& bumpShaders):
        Material(BSDFType(BSDFGlossy | BSDFReflection), bumpShaders), 
        mGlossyFactor(Kg), mExp(exponent), mEta(index),
        mK(absorption), mFresnelType(Conductor) {}


    class TransparentMaterial : public Material {
    public:
        TransparentMaterial(const ColorTexturePtr& Kr, const ColorTexturePtr& Kt, 
            float index, const BumpShaders& bumpShaders = BumpShaders());
        Color bsdf(const Fragment& fragment, const Vector3& wo,
            const Vector3& wi, BSDFType type) const;

        // given wo sample a wi and corresponding bsdf value to return
        Color sampleBSDF(const Fragment& fragment, 
            const Vector3& wo, const BSDFSample& bsdfSample, Vector3* wi, 
            float* pdf, BSDFType type, BSDFType* sampledType) const;

        float pdf(const Fragment& fragment, 
            const Vector3& wo, const Vector3& wi, BSDFType type) const;

    private:
        ColorTexturePtr mReflectFactor;
        ColorTexturePtr mRefractFactor;
        float mEtai;
        float mEtat;
    };

    inline TransparentMaterial::TransparentMaterial(const ColorTexturePtr& Kr,
        const ColorTexturePtr& Kt, float index, const BumpShaders& bumpShaders):
        Material(BSDFType(BSDFSpecular | BSDFReflection | BSDFTransmission),
        bumpShaders), mReflectFactor(Kr), mRefractFactor(Kt), 
        mEtai(1.0f), mEtat(index) {}

    // there is only one possible wi for specified wo, specular reflection
    // is a delta distribution function, we count on sample way to get
    // the wi by feeding in wo
    inline Color TransparentMaterial::bsdf(const Fragment& fragment, 
        const Vector3& wo, const Vector3& wi, BSDFType type) const {
        return Color::Black;
    }

    inline float TransparentMaterial::pdf(const Fragment& fragment,
        const Vector3& wo, const Vector3& wi, BSDFType type) const {
        return 0.0f;
    }


    class MirrorMaterial : public Material {
    public:
        MirrorMaterial(const ColorTexturePtr& Kr, float index, 
            float absorption, const BumpShaders& bumpShaders = BumpShaders());
        Color bsdf(const Fragment& fragment, const Vector3& wo,
            const Vector3& wi, BSDFType type) const;

        // given wo sample a wi and corresponding bsdf value to return
        Color sampleBSDF(const Fragment& fragment, 
            const Vector3& wo, const BSDFSample& bsdfSample, Vector3* wi, 
            float* pdf, BSDFType type, BSDFType* sampledType) const;

        float pdf(const Fragment& fragment, 
            const Vector3& wo, const Vector3& wi, BSDFType type) const;

    private:
        ColorTexturePtr mReflectFactor;
        float mEta;
        float mK;
    };

    inline MirrorMaterial::MirrorMaterial(const ColorTexturePtr& Kr, 
        float index, float absorption, const BumpShaders& bumpShaders):
        Material(BSDFType(BSDFSpecular | BSDFReflection), bumpShaders), 
        mReflectFactor(Kr), mEta(index), mK(absorption) {}

    inline Color MirrorMaterial::bsdf(const Fragment& fragment, 
        const Vector3& wo, const Vector3& wi, BSDFType type) const {
        return Color::Black;
    }

    inline float MirrorMaterial::pdf(const Fragment& fragment,
        const Vector3& wo, const Vector3& wi, BSDFType type) const {
        return 0.0f;
    }


    class SubsurfaceMaterial : public Material {
    public:
        SubsurfaceMaterial(const ColorTexturePtr& absorb, 
            const ColorTexturePtr& scatterPrime, 
            const ColorTexturePtr& Kr, float eta, float g,
            const BumpShaders& bumpShaders = BumpShaders());

        SubsurfaceMaterial(const Color& Kd, 
            const Color& diffuseMeanFreePath, 
            const ColorTexturePtr& Kr, 
            float eta, float g,
            const BumpShaders& bumpShaders = BumpShaders());

        ~SubsurfaceMaterial();

        Color bsdf(const Fragment& fragment, const Vector3& wo,
            const Vector3& wi, BSDFType type) const;

        Color sampleBSDF(const Fragment& fragment, 
            const Vector3& wo, const BSDFSample& bsdfSample, Vector3* wi, 
            float* pdf, BSDFType type, BSDFType* sampledType) const;

        float pdf(const Fragment& fragment, 
            const Vector3& wo, const Vector3& wi, BSDFType type) const;

        const BSSRDF* getBSSRDF() const;

    private:
        ColorTexturePtr mReflectFactor;
        float mEta;
        BSSRDF* mBSSRDF;
    };

    inline SubsurfaceMaterial::SubsurfaceMaterial(
        const ColorTexturePtr& absorb, const ColorTexturePtr& scatterPrime, 
        const ColorTexturePtr& Kr, float eta, float g, 
        const BumpShaders& bumpShaders): 
        Material(BSDFAll, bumpShaders), mReflectFactor(Kr), mEta(eta) {
        mBSSRDF = new BSSRDF(absorb, scatterPrime, eta, g);
    }

    inline SubsurfaceMaterial::SubsurfaceMaterial(const Color& Kd, 
        const Color& diffuseMeanFreePath, const ColorTexturePtr& Kr, 
        float eta, float g, const BumpShaders& bumpShaders):
        Material(BSDFAll, bumpShaders), mReflectFactor(Kr), mEta(eta) {
        mBSSRDF = new BSSRDF(Kd, diffuseMeanFreePath, eta, g);
    }

    inline SubsurfaceMaterial::~SubsurfaceMaterial() {
        if(mBSSRDF != NULL) {
            delete mBSSRDF;
            mBSSRDF = NULL;
        }
    }

    inline Color SubsurfaceMaterial::bsdf(const Fragment& fragment, 
        const Vector3& wo, const Vector3& wi, BSDFType type) const {
        return Color::Black;
    }

    inline float SubsurfaceMaterial::pdf(const Fragment& fragment,
        const Vector3& wo, const Vector3& wi, BSDFType type) const {
        return 0.0f;
    }

    inline const BSSRDF* SubsurfaceMaterial::getBSSRDF() const {
        return mBSSRDF;
    }


    class MaskMaterial : public Material {
    public:
        MaskMaterial(const FloatTexturePtr& alphaMask, 
            const ColorTexturePtr& transparentColor,
            const MaterialPtr& maskedMaterial);

        Color bsdf(const Fragment& fragment, const Vector3& wo,
            const Vector3& wi, BSDFType type) const;

        Color sampleBSDF(const Fragment& fragment, 
            const Vector3& wo, const BSDFSample& bsdfSample, Vector3* wi, 
            float* pdf, BSDFType type, BSDFType* sampledType) const;

        float pdf(const Fragment& fragment, 
            const Vector3& wo, const Vector3& wi, BSDFType type) const;

        // override the bump mapping since it's the masked material
        // should do the job
        void perturb(Fragment* fragment) const;
    private:
        FloatTexturePtr mAlphaMask;
        ColorTexturePtr mTransparentColor;
        MaterialPtr mMaskedMaterial;
    };

    inline MaskMaterial::MaskMaterial(const FloatTexturePtr& alphaMask, 
        const ColorTexturePtr& transparentColor,
        const MaterialPtr& maskedMaterial): 
        Material(BSDFType(maskedMaterial->getType() | BSDFNull), 
        BumpShaders()), mAlphaMask(alphaMask), 
        mTransparentColor(transparentColor),
        mMaskedMaterial(maskedMaterial) {}

    inline void MaskMaterial::perturb(Fragment* fragment) const {
        mMaskedMaterial->perturb(fragment);
    }


    class LambertMaterialCreator : public 
        Creator<Material , const ParamSet&, const SceneCache&> {
    public:
        Material* create(const ParamSet& params, 
            const SceneCache& sceneCache) const;
    };

    
    class BlinnMaterialCreator : public 
        Creator<Material , const ParamSet&, const SceneCache&> {
    public:
        Material* create(const ParamSet& params, 
            const SceneCache& sceneCache) const;
    };


    class TransparentMaterialCreator : public 
        Creator<Material , const ParamSet&, const SceneCache&> {
    public:
        Material* create(const ParamSet& params, 
            const SceneCache& sceneCache) const;
    };


    class MirrorMaterialCreator : public 
        Creator<Material , const ParamSet&, const SceneCache&> {
    public:
        Material* create(const ParamSet& params, 
            const SceneCache& sceneCache) const;
    };


    class SubsurfaceMaterialCreator : public 
        Creator<Material , const ParamSet&, const SceneCache&> {
    public:
        Material* create(const ParamSet& params, 
            const SceneCache& sceneCache) const;
    };


    class MaskMaterialCreator : public 
        Creator<Material , const ParamSet&, const SceneCache&> {
    public:
        Material* create(const ParamSet& params, 
            const SceneCache& sceneCache) const;
    };
}

#endif //GOBLIN_MATERIAL_H
