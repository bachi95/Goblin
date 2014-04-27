#ifndef GOBLIN_MATERIAL_H
#define GOBLIN_MATERIAL_H

#include "GoblinColor.h"
#include "GoblinTexture.h"
#include "GoblinUtils.h"

namespace Goblin {
    struct Vertex;
    class Fragment;
    class Vector3;
    class Matrix3;
    class Sampler;
    class Sample;
    
    enum BSDFType {
        BSDFReflection = 1 << 0,
        BSDFTransmission = 1 << 1,
        BSDFDiffuse = 1 << 2,
        BSDFGlossy = 1 << 3,
        BSDFSpecular = 1 << 4,
        BSDFAll = BSDFReflection | BSDFTransmission |
            BSDFDiffuse | BSDFGlossy | BSDFSpecular
    };

    struct BSDFSampleIndex {
        BSDFSampleIndex() {}
        BSDFSampleIndex(Sampler* sampler, int requestNum);
        uint32_t samplesNum;
        uint32_t directionIndex;
    };

    struct BSDFSample {
        BSDFSample();
        BSDFSample(const Sample& sample,
            const BSDFSampleIndex& index, uint32_t n);
        float uDirection[2];
    };

    // This serves the purpose of what GPU shader usually do
    // perturb the geometry info(position, normal, tangent....)
    // calculate the bsdf value, hold the texture reference...etc
    class Material {
    public:
        virtual ~Material() {}

        virtual Color bsdf(const Fragment& fragment, const Vector3& wo, 
            const Vector3& wi, BSDFType type = BSDFAll) const = 0;

        virtual Color sampleBSDF(const Fragment& fragment, 
            const Vector3& wo, const BSDFSample& bsdfSample, Vector3* wi, 
            float* pdf, BSDFType type = BSDFAll, 
            BSDFType* sampledType = NULL) const = 0;

        virtual float pdf(const Fragment& fragment,
            const Vector3& wo, const Vector3& wi,
            BSDFType type = BSDFAll) const = 0;

    protected:
        bool matchType(BSDFType type, BSDFType toMatch) const;

        BSDFType getType(const Vector3& wo, const Vector3& wi, 
            const Fragment& f, BSDFType type) const;

        bool sameHemisphere(const Fragment& fragment,
            const Vector3& wo, const Vector3& wi) const;

        float specularReflect(const Fragment& fragment, const Vector3& wo,
            Vector3* wi, float etai, float etat) const;

        float specularRefract(const Fragment& fragment, const Vector3& wo,
            Vector3* wi, float etai, float etat) const;

        // material util to get the fresnel factor
        float fresnelDieletric(float cosi, float etai, float etat) const;
    };

    inline bool Material::matchType(BSDFType type, BSDFType toMatch) const {
        return (type & toMatch) == toMatch;
    }

    typedef boost::shared_ptr<Material> MaterialPtr;

    class LambertMaterial : public Material {
    public:
        LambertMaterial(const TexturePtr& Kd);
        Color bsdf(const Fragment& fragment, const Vector3& wo, 
            const Vector3& wi, BSDFType type) const;

        Color sampleBSDF(const Fragment& fragment, 
            const Vector3& wo, const BSDFSample& bsdfSample, Vector3* wi, 
            float* pdf, BSDFType type, BSDFType* sampledType) const;

        float pdf(const Fragment& fragment, 
            const Vector3& wo, const Vector3& wi, BSDFType type) const;

    private:
        TexturePtr mDiffuseFactor;
    };

    inline LambertMaterial::LambertMaterial(const TexturePtr& Kd):
        mDiffuseFactor(Kd) {}


    class TransparentMaterial : public Material {
    public:
        TransparentMaterial(const TexturePtr& Kr, const TexturePtr& Kt, 
            float index);
        Color bsdf(const Fragment& fragment, const Vector3& wo,
            const Vector3& wi, BSDFType type) const;

        // given wo sample a wi and corresponding bsdf value to return
        Color sampleBSDF(const Fragment& fragment, 
            const Vector3& wo, const BSDFSample& bsdfSample, Vector3* wi, 
            float* pdf, BSDFType type, BSDFType* sampledType) const;

        float pdf(const Fragment& fragment, 
            const Vector3& wo, const Vector3& wi, BSDFType type) const;

    private:
        TexturePtr mReflectFactor;
        TexturePtr mRefractFactor;
        float mEtai;
        float mEtat;
    };

    inline TransparentMaterial::TransparentMaterial(const TexturePtr& Kr,
        const TexturePtr& Kt, float index):
        mReflectFactor(Kr), mRefractFactor(Kt), mEtai(1.0f), mEtat(index) {}

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
}

#endif //GOBLIN_MATERIAL_H
