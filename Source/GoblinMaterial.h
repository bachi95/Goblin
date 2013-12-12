#ifndef GOBLIN_MATERIAL_H
#define GOBLIN_MATERIAL_H

#include "GoblinUtils.h"
#include "GoblinColor.h"

namespace Goblin {
    struct Vertex;
    struct Fragment;
    class Vector3;
    class Matrix3;
    
    enum BSDFType {
        BSDFReflection = 1 << 0,
        BSDFTransmission = 1 << 1,
        BSDFDiffuse = 1 << 2,
        BSDFGlossy = 1 << 3,
        BSDFSpecular = 1 << 4,
        BSDFAll = BSDFReflection | BSDFTransmission |
            BSDFDiffuse | BSDFGlossy | BSDFSpecular
    };

    // This serves the purpose of what GPU shader usually do
    // perturb the geometry info(position, normal, tangent....)
    // calculate the bsdf value, hold the texture reference...etc
    class Material {
    public:
        virtual ~Material() {}

        virtual Color bsdf(const Fragment& fragment, const Vector3& wo, 
            const Vector3& wi, BSDFType type = BSDFAll) const = 0;

        virtual Color sampleBSDF(const Fragment& fragment, const Vector3& wo,
            Vector3* wi, BSDFType type = BSDFAll) const = 0;

    protected:
        bool matchType(BSDFType type, BSDFType toMatch) const;

        BSDFType getType(const Vector3& wo, const Vector3& wi, 
            const Fragment& f, BSDFType type) const;

        // get the matrix convert world space vector to shading space
        // formed by tangent, bitangent, normal
        // normal is Vector3(0, 0, 1) in this shading space
        Matrix3 getShadingMatrix(const Fragment& fragment) const;

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
        LambertMaterial(const Color& Kd);
        Color bsdf(const Fragment& fragment, const Vector3& wo, 
            const Vector3& wi, BSDFType type) const;

        Color sampleBSDF(const Fragment& fragment, const Vector3& wo,
            Vector3* wi, BSDFType type = BSDFAll) const;

    private:
        Color mDiffuseFactor;
    };

    inline LambertMaterial::LambertMaterial(const Color& Kd):
        mDiffuseFactor(Kd) {}


    class TransparentMaterial : public Material {
    public:
        TransparentMaterial(const Color& Kr, const Color& Kt, float index);
        Color bsdf(const Fragment& fragment, const Vector3& wo,
            const Vector3& wi, BSDFType type) const;

        // given wo sample a wi and corresponding bsdf value to return
        Color sampleBSDF(const Fragment& fragment, const Vector3& wo, 
            Vector3* wi, BSDFType type = BSDFAll) const;

    private:
        Color mReflectFactor;
        Color mRefractFactor;
        float mEtai;
        float mEtat;
    };

    inline TransparentMaterial::TransparentMaterial(const Color& Kr,
        const Color& Kt, float index):
        mReflectFactor(Kr), mRefractFactor(Kt), mEtai(1.0f), mEtat(index) {}

    // there is only one possible wi for specified wo, specular reflection
    // is a delta distribution function, we are count on sample way to get
    // the wi by feeding in wo
    inline Color TransparentMaterial::bsdf(const Fragment& fragment, 
        const Vector3& wo, const Vector3& wi, BSDFType type) const {
        return Color::Black;
    }
}

#endif //GOBLIN_MATERIAL_H
