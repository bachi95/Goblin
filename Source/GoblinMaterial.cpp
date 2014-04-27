#include "GoblinMaterial.h"
#include "GoblinVector.h"
#include "GoblinMatrix.h"
#include "GoblinVertex.h"
#include "GoblinGeometry.h"
#include "GoblinSampler.h"

namespace Goblin {

    BSDFSampleIndex::BSDFSampleIndex(Sampler* sampler, 
        int requestNum) {
        SampleIndex twoDIndex = sampler->requestTwoDQuota(requestNum);
        samplesNum = twoDIndex.sampleNum;
        directionIndex = twoDIndex.offset;
    }

    BSDFSample::BSDFSample() {
        uDirection[0] = randomFloat();
        uDirection[1] = randomFloat();
    }

    BSDFSample::BSDFSample(const Sample& sample,
        const BSDFSampleIndex& index, uint32_t n) {
        uDirection[0] = sample.u2D[index.directionIndex][2 * n];
        uDirection[1] = sample.u2D[index.directionIndex][2 * n + 1];
    }

    BSDFType Material::getType(const Vector3& wo, const Vector3& wi,
        const Fragment& fragment, BSDFType type) const {
        const Vector3& n = fragment.getNormal();
        if(dot(n, wo) * dot(n, wi) > 0.0f) {
            type = BSDFType(type & ~BSDFTransmission);
        } else {
            type = BSDFType(type & ~BSDFReflection);
        }
        return type;
    }

    bool Material::sameHemisphere(const Fragment& fragment,
        const Vector3& wo, const Vector3& wi) const {
        const Vector3& normal = fragment.getNormal();
        return dot(wo, normal) * dot(wi, normal) > 0.0f;
    }

    float Material::specularReflect(const Fragment& fragment, 
        const Vector3& wo, Vector3* wi, float etai, float etat) const {
        Vector3 n = fragment.getNormal();
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
        Vector3 n = fragment.getNormal();
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
            f += mDiffuseFactor->lookup(fragment) * INV_PI;
        }
        return f;
    }

    Color LambertMaterial::sampleBSDF(const Fragment& fragment, 
        const Vector3& wo, const BSDFSample& bsdfSample, Vector3* wi, 
        float* pdf, BSDFType type, BSDFType* sampledType) const {
        BSDFType materialType = BSDFType(BSDFDiffuse | BSDFReflection);
        if(!matchType(type, materialType)) {
            *pdf = 0.0f;
            return Color::Black;
        }
        float u1 = bsdfSample.uDirection[0];
        float u2 = bsdfSample.uDirection[1];
        Vector3 wiLocal = cosineSampleHemisphere(u1, u2);
        // flip it if wo and normal at different sides of hemisphere
        if(dot(wo, fragment.getNormal()) < 0.0f) {
            wiLocal *= -1.0f;
        }
        Matrix3 shadeToWorld = fragment.getWorldToShade().transpose();
        *wi = shadeToWorld * wiLocal;;
        *pdf = this->pdf(fragment, wo, *wi, materialType);
        if(sampledType) {
            *sampledType = materialType;
        }
        return mDiffuseFactor->lookup(fragment) * INV_PI;
    }

    float LambertMaterial::pdf(const Fragment& fragment,
        const Vector3& wo, const Vector3& wi,
        BSDFType type) const {
        if(!matchType(type, BSDFType(BSDFDiffuse | BSDFReflection))) {
            return 0.0f;
        }
        return sameHemisphere(fragment, wo, wi)? 
            absdot(fragment.getNormal(), wi) * INV_PI : 0.0f;
    }

    Color TransparentMaterial::sampleBSDF(const Fragment& fragment, 
        const Vector3& wo, const BSDFSample& bsdfSample, Vector3* wi, 
        float* pdf, BSDFType type, BSDFType* sampledType) const {
        int nMatch = 0;
        if(matchType(type, BSDFType(BSDFSpecular | BSDFReflection))) {
            nMatch++;
        }
        if(matchType(type, BSDFType(BSDFSpecular | BSDFTransmission))) {
            nMatch++;
        }

        Color f(Color::Black);
        if(nMatch == 1) {
            if(matchType(type, BSDFReflection)) {
                f = mReflectFactor->lookup(fragment) * 
                    specularReflect(fragment, wo, wi, mEtai, mEtat);
                if(sampledType) {
                    *sampledType = BSDFType(BSDFSpecular | BSDFReflection);
                }
            } else {
                f = mRefractFactor->lookup(fragment) *
                    specularRefract(fragment, wo, wi, mEtai, mEtat);
                if(sampledType) {
                    *sampledType = BSDFType(BSDFSpecular | BSDFTransmission);
                }
            }
            *pdf = 1.0f;
        } else if(nMatch == 2) {
            Vector3 wReflect;
            Vector3 wRefract;
            float reflect = specularReflect(fragment, wo, 
                &wReflect, mEtai, mEtat);
            float refract = specularRefract(fragment, wo, 
                &wRefract, mEtai, mEtat);
            // use fresnel factor to do importance sampling
            float fresnel = reflect * absdot(wReflect, fragment.getNormal());
            float reflectChance = fresnel;
            bool doReflect = randomFloat() < reflectChance;

            if(doReflect) {
                f = mReflectFactor->lookup(fragment) * reflect;
                *wi = wReflect;
                if(sampledType) {
                    *sampledType = BSDFType(BSDFSpecular | BSDFReflection);
                }
                *pdf = reflectChance;
            } else {
                f = mRefractFactor->lookup(fragment) * refract;
                    specularRefract(fragment, wo, wi, mEtai, mEtat);
                *wi = wRefract;
                if(sampledType) {
                    *sampledType = BSDFType(BSDFSpecular | BSDFTransmission);
                }
                *pdf = 1.0f - reflectChance;
            }
        } else {
            *pdf = 0.0f;
        }
        return f;
    }

}
