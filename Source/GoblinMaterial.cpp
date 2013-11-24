#include "GoblinMaterial.h"
#include "GoblinVector.h"
#include "GoblinVertex.h"

namespace Goblin {
    Color LambertMaterial::bsdf(const Vertex& sg, const Vector3& wi, 
        const Vector3& wo, BSDFType type) const {
        Color f(Color::Black);
        if(type & BSDFReflection) {
            f += mDiffuseFactor * INV_PI;
        }
        return f;
    }
}