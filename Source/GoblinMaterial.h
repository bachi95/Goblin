#ifndef GOBLIN_MATERIAL_H
#define GOBLIN_MATERIAL_H

#include "GoblinUtils.h"
#include "GoblinColor.h"

namespace Goblin {
    struct Vertex;
    class Vector3;

    enum BSDFType {
        BSDFReflection = 1 << 0,
        BSDFTransmission = 1 << 1,
        BSDFAll = BSDFReflection | BSDFTransmission
    };
    // This serves the purpose of what GPU shader usually do
    // perturb the geometry info(position, normal, tangent....)
    // calculate the bsdf value, hold the texture reference...etc
    class Material {
    public:
        virtual ~Material() {}
        virtual Color bsdf(const Vertex& sg, const Vector3& wi, 
            const Vector3& wo, BSDFType type = BSDFAll) const = 0;
    };

    typedef boost::shared_ptr<Material> MaterialPtr;

    class LambertMaterial : public Material {
    public:
        LambertMaterial(const Color& Kd);
        Color bsdf(const Vertex& sg, const Vector3& wi, 
            const Vector3& wo, BSDFType type) const;
    private:
        Color mDiffuseFactor;
    };

    inline LambertMaterial::LambertMaterial(const Color& Kd):
        mDiffuseFactor(Kd) {}
}

#endif //GOBLIN_MATERIAL_H
