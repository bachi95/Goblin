#include "GoblinMaterial.h"
#include "GoblinVector.h"
#include "GoblinMatrix.h"
#include "GoblinVertex.h"
#include "GoblinGeometry.h"

namespace Goblin {

    BSDFType Material::getType(const Vector3& wo, const Vector3& wi,
        const Fragment& fragment, BSDFType type) const {
        const Vector3& n = fragment.normal;
        if(dot(n, wo) * dot(n, wi) > 0.0f) {
            type = BSDFType(type & ~BSDFTransmission);
        } else {
            type = BSDFType(type & ~BSDFReflection);
        }
        return type;
    }

    Matrix3 Material::getShadingMatrix(const Fragment& fragment) const {
        Vector3 n = fragment.normal;
        // dpdu is not necessarily orthgnal to normal, use gram-shmidt
        // to tick out normal component from dpdu to form tangent
        Vector3 t = normalize(fragment.dpdu - n * dot(fragment.dpdu, n));
        Vector3 b = cross(n, t);
        return Matrix3(
            t.x, t.y, t.z,
            b.x, b.y, b.z,
            n.x, n.y, n.z);
    }

    float Material::specularReflect(const Fragment& fragment, 
        const Vector3& wo, Vector3* wi, float etai, float etat) const {
        Vector3 n = fragment.normal;
        // wo(face outward) is the input ray(with n form angle i), 
        // wi is the reflect ray
        float cosi = dot(n, wo);
        float ei = etai;
        float et = etat;
        bool entering = cosi > 0.0f;
        if(!entering) {
            swap(ei, et);
            n = -n;
            cosi = -cosi;
        }
        float f = fresnelDieletric(cosi, ei, et);
        // Wi = Wo - 2 * WoPerpN = Wo - 2(Wo - dot(N, Wo))N =
        // 2(dot(N, Wo))N - Wo
        *wi = 2 * cosi * n - wo;
        // specular reflection angle should be the same as input angle
        float cosr = cosi;
        return f / cosr;
    }

    float Material::specularRefract(const Fragment& fragment,
        const Vector3& wo, Vector3* wi, float etai, float etat) const {
        // wo(face outward) is the input ray(with n form angle i), 
        // wi is the refract ray(with -n form angle t)
        Vector3 n = fragment.normal;
        // wo(face outward) is the input ray(with n form angle i), 
        // wi is the reflect ray
        float cosi = dot(n, wo);
        float ei = etai;
        float et = etat;
        bool entering = cosi > 0.0f;
        if(!entering) {
            swap(ei, et);
            n = -n;
            cosi = -cosi;
        }
        float f = fresnelDieletric(cosi, ei, et);
        // total reflection
        if(f == 1.0f) {
            return 0.0f;
        }
        //Wi = -N * cost - WoPerpN * sint / sini =
        //-N * cost - (sint / sini) * (Wo - dot(N, Wo) * N) =
        //-N * sqrt(1 -(etai/etat)^2 * (1 - (dot(N, Wo))^2))) +
        //etai/etat(Wo - dot(N, Wo) * N) =
        //N * (etai/etat * dot(N, Wo) - 
        //sqrt(1 - (etai/etat)^2(1 - (dot(N, Wo))^2)) -
        //etai / etat * Wo
        float eta = ei / et;
        *wi = normalize(n * (eta * cosi - 
            sqrt(max(0.0f, 1.0f - eta * eta * (1.0f - cosi * cosi)))) - 
            eta * wo);
        return (1.0f - f) / absdot(*wi, n);
    }

    // close approximation fresnel for dieletric material
    float Material::fresnelDieletric(float cosi, float etai, 
        float etat) const {
        cosi = clamp(cosi, -1.0f, 1.0f);
        float sint = (etai / etat) * sqrt(max(0.0f, 1.0f - cosi * cosi));
        // total reflection
        if(sint >= 1.0f) {
            return 1.0f;
        }
        float cost = sqrt(max(0.0f, 1 - sint * sint));
        cosi = fabs(cosi);

        float rParl = ((etat * cosi) - (etai * cost)) /
            ((etat * cosi) + (etai * cost));
        float rPerp = ((etai * cosi) - (etat * cost)) /
            ((etai * cosi) + (etat * cost));
        return (rParl * rParl + rPerp * rPerp) / 2.0f;
    }

    Color LambertMaterial::bsdf(const Fragment& fragment, const Vector3& wo, 
        const Vector3& wi, BSDFType type) const {
        Color f(Color::Black);
        type = getType(wo, wi, fragment, type);
        if(matchType(type, BSDFType(BSDFDiffuse | BSDFReflection))) {
            f += mDiffuseFactor * INV_PI;
        }
        return f;
    }

    Color LambertMaterial::sampleBSDF(const Fragment& fragment, 
        const Vector3& wo, Vector3* wi, BSDFType type) const {
        Color f(Color::Black);
        if(matchType(type, BSDFType(BSDFDiffuse | BSDFReflection))) {
            f += mDiffuseFactor * INV_PI;
        }
        //TODO wi need to be determined by monte carlo sampling
        return f;
    }

    Color TransparentMaterial::sampleBSDF(const Fragment& fragment, 
        const Vector3& wo, Vector3* wi, BSDFType type) const {
        if(matchType(type, BSDFType(BSDFSpecular | BSDFReflection))) {
            Color f = mReflectFactor * 
                specularReflect(fragment, wo, wi, mEtai, mEtat);
            return f;
        }
        if(matchType(type, BSDFType(BSDFSpecular | BSDFTransmission))) {
            Color f = mRefractFactor * 
                specularRefract(fragment, wo, wi, mEtai, mEtat);
            return f;
        }
        return Color::Black;
    }

}
