#include "GoblinMaterial.h"
#include "GoblinVector.h"
#include "GoblinMatrix.h"
#include "GoblinVertex.h"
#include "GoblinGeometry.h"
#include "GoblinSampler.h"

namespace Goblin {

    BSDFSampleIndex::BSDFSampleIndex(SampleQuota* sampleQuota, 
        int requestNum) {
        SampleIndex twoDIndex = sampleQuota->requestTwoDQuota(requestNum);
        samplesNum = twoDIndex.sampleNum;
        directionIndex = twoDIndex.offset;
    }

    BSDFSample::BSDFSample(const RNG& rng) {
        uDirection[0] = rng.randomFloat();
        uDirection[1] = rng.randomFloat();
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

    void Material::perturb(Fragment* fragment) const {
        /*
         * let bumpmap(u, v) be D
         * p'(u, v) = p(u, v) + n(u, v) * D
         * we want to get the new normal n'(u, v),
         * which is cross(dp'du, dp'dv)
         * dp'du = d(p'(u, v)/du = dp/du + dD/du * n d+ dn/du * D
         * we have dp/du already
         * dn/du * D can be ignored when D is small
         * dD/du can be approximated by forward differential approximation
         * (D(u + du, v) - D(u, v)) / du
         * dp'dv is calculated with the same way above
         */
        if(mBumpMap) {
            Vector3 p = fragment->getPosition();
            Vector3 n = fragment->getNormal();
            Vector2 uv = fragment->getUV();
            float bumpD = mBumpMap->lookup(*fragment);

            float du = 0.002f;
            Fragment fdu = *fragment; 
            fdu.setPosition(p + du * fragment->getDPDU());
            fdu.setUV(uv + Vector2(du, 0.0f));
            float bumpDdu = mBumpMap->lookup(fdu);
            Vector3 bumpDPDU = fragment->getDPDU() + 
                (bumpDdu - bumpD) / du * n;

            float dv = 0.002f;
            Fragment fdv = *fragment;
            fdv.setPosition(p + dv * fragment->getDPDV());
            fdv.setUV(uv + Vector2(0.0f, dv));
            float bumpDdv = mBumpMap->lookup(fdv);
            Vector3 bumpDPDV = fragment->getDPDV() +
                (bumpDdv - bumpD) / dv * n;

            Vector3 bumpN = normalize(cross(bumpDPDU, bumpDPDV));
            // case that intersect back face
            if(dot(bumpN, n) < 0.0f) {
                bumpN *= -1.0f;
            }
            // now set back the bumped value
            fragment->setNormal(bumpN);
            fragment->setDPDU(bumpDPDU);
            fragment->setDPDV(bumpDPDV);
        }
        return;
    }

    bool Material::sameHemisphere(const Fragment& fragment,
        const Vector3& wo, const Vector3& wi) const {
        const Vector3& normal = fragment.getNormal();
        return dot(wo, normal) * dot(wi, normal) > 0.0f;
    }

    float Material::specularReflectDieletric(const Fragment& fragment, 
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

    float Material::specularReflectConductor(const Fragment& fragment, 
        const Vector3& wo, Vector3* wi, float eta, float k) const {
        Vector3 n = fragment.getNormal();
        float cosi = dot(n, wo);
        // back facing
        if(cosi <= 0.0f) {
            return 0.0f;
        }
        float f = fresnelConductor(cosi, eta, k);
        *wi = 2 * cosi * n - wo;
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
        /*
         * Wi = -N * cost - WoPerpN * sint / sini =
         * -N * cost - (sint / sini) * (Wo - dot(N, Wo) * N) =
         * -N * sqrt(1 - (etai / etat)^2 * (1 - (dot(N, Wo))^2))) +
         * etai / etat * (Wo - dot(N, Wo) * N) =
         * N * (etai/etat * dot(N, Wo) - 
         * sqrt(1 - (etai/etat)^2(1 - (dot(N, Wo))^2)) -
         * etai / etat * Wo
         */
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

    float Material::fresnelConductor(float cosi, float eta, float k) const {
        float tmp = (eta * eta + k * k);
        float cosi2 = cosi * cosi;
        float rParl2 = (tmp * cosi2 - 2.0f * eta * cosi + 1.0f) /
            (tmp * cosi2 + 2.0f * eta * cosi + 1.0f);
        float rPerp2 = (tmp - 2.0f * eta * cosi + cosi2) /
            (tmp + 2.0f * eta * cosi + cosi2);
        return (rParl2 + rPerp2) * 0.5f;
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

    /*
     * implementation based on Torrance-Sparrow microfacet model with
     * Blinn microfacet distribution
     * Let wh = wi + wo (half-angle vector)
     * Define D(wh) the probablility distribution function that gives
     * possiblity microfacet facing wh
     * 
     * Differential Flux received by microfacets with half angle wh:
     * dFlux = Li(wi) * dwi  * cos(thetah) * dA(wh)
     * where dA(wh) is the area measure of the microfacets with orientation wh
     * dA(wh) = (D(wh) * dwh) * dA =>
     * dFlux = Li(wi) * dwi * cos(thetah) * D(wh) * dwh * dA
     * 
     * since in Torrance-Sparrow model we assume each microfacet 
     * perfect specular, we apply fresnel factor on dFlux to get outgoing Flux
     * dFluxO = Fr(wo)dFlux =>
     * L(wo) = dFluxO / (dwo * cos(thetao) * dA) =
     *     Fr(wo) * Li(wi) * dwi  * D(wh) * dwh * dA * cos(thetah) / 
     *     (dwo * cos(thetao) * dA)
     * 
     * dwh / dwo = (sin(thetah) * dthetah * dphih) / 
     *     (sin(thetao) * dthetao * dphio)
     * because wo is computed by reflect wi about wh, thetao = 2thetah =>
     * dwh / dwo = (sin(thetah) * dthetah * dphih) /
     *     (sin(2thetah) * 2dthetah * dphio) =
     *     1 / 4cos(thetah) = 1 / (4 * dot(wo, wh)) = 1 / (4 * dot(wi, wh)) =>
     * dwh = dwo / 4cos(thetah)
     * we can then substitute the above L(wo) to get: 
     * L(wo) = Fr(wo) * Li(wi) * D(wh) * dwi / (4 * cos(thetao)) =>
     * bsdf(wo, wi) = D(wh) * Fr(wo) / (4 * cos(thetao) * cos(thetai))
     * 
     * There is another geometric attenuation term that describes the fraction
     * of reflected light without maked or shadowed by neighbor microfacet,
     * the derivation assume microfacets as infinitely long v-shaped grooves, 
     * the detail derivation can be seen at:
     * Blinn J. F. 1977. 
     * Models of light reflection for computer synthesized pictures
     * we'll denote this geometry attenuation term as
     * G(wo, wi) = min(
     *     1, 
     *     2 * dot(n, wh) * dot(n, wo) / dot(wo, wh),
     *     2 * dot(n, wh) * dot(n, wi) / dot(wo, wh) )
     *
     * bsdf(wo, wi) = D(wh) * G(wo, wi) * Fr(wo) / 
     *     (4 * cos(thetao) * cos(thetai))
     *
     * The D(wh) is a exponential falloff distribution
     * D(wh) = c * pow(dot(wh, n), e)
     * To make this distribution energy conserved, the integration of D(wh)
     * over dwh projected to dA should be dA:
     * integragte(c * pow(dot(wh, n), e) * cos(thetah), dwh) over hemisphere = 1
     * 2 * c * PI * integrate(pow(cos(thetah), e + 1) * sin(thetah), dtheta) 
     * over 0 to PI / 2 = 1 => substute cos(thetah) with u:
     * 2 * c * PI * integrate(pow(u, e + 1), du) over 0 to 1 = 1 =>
     * 2 * c * PI / (e + 2) = 1 => c = (e + 2) / 2PI =>
     * D(wh) = (e + 2) / 2PI * pow(dot(wh, n), e)
     *
     */
    Color BlinnMaterial::bsdf(const Fragment& fragment, const Vector3& wo, 
        const Vector3& wi, BSDFType type) const {
        type = getType(wo, wi, fragment, type);
        if(matchType(type, BSDFType(BSDFGlossy | BSDFReflection))) {
            Vector3 n = fragment.getNormal();
            float cosi = absdot(n, wi);
            float coso = absdot(n, wo);
            if(cosi == 0.0f || coso == 0.0f) {
                return Color::Black;
            }
            Vector3 wh = normalize(wo + wi);
            float cosh = absdot(n, wh);
            float exp = mExp->lookup(fragment);
            // blinn distribution
            float D = (exp + 2.0f) * INV_TWOPI * 
                pow(cosh, exp);
            float woDotWh = absdot(wo, wh);
            // geometry attenuation term
            float G = min(1.0f, min(2.0f * cosh * coso / woDotWh, 
                2.0f * cosh * cosi / woDotWh));
            // fresnel factor
            float F = 1.0f;
            if(mFresnelType == Dieletric) {
                F = fresnelDieletric(woDotWh, 1.0f, mEta);
            } else if(mFresnelType == Conductor) {
                F = fresnelConductor(woDotWh, mEta, mK);
            }
            return mGlossyFactor->lookup(fragment) * D * G * F / 
                (4.0f * cosi * coso);
        }
        return Color::Black;
    }

    /*
     * first sample wh with Blinn distribution, then compute wi by reflect
     * wo over wh
     * pdf(thetah, phih) = c * D(wh) = c * pow(cos(thetah), e)
     * integrate(D(wh), dw) over hemisphere = 1 =>
     * 2 * c * PI * integrate(pow(cos(thetah), e) * sin(thetah), dw) over 
     * 0-PI / 2 = 1 => substitute cos(thetah) with u
     * 2 * c * PI * integrate(pow(u, e), du) over 0-1 = 1 =>
     * c = (e + 1) / 2PI =>
     * pdf(thetah, phih) = (e + 1) * pow(cos(thetah), e) / 2PI
     * since phih doesn't affect D, pdf(phih) is a separable uniform pdf
     * pdf(phih) = 1 / 2PI
     * pdf(thetah) = pdf(cos(thetah)) = (e + 1) * pow(cos(tehtah), e)
     * cdf(cos(thetah)) = pow(cos(tehtah), e + 1)
     * inverse method:
     * cos(thetah) = pow(u1, 1 / (e + 1))
     * phi = 2PI * u2
     * with the above cos(thetah) and phi we can get wh in shading space
     * we can then transform it back to world space and get wi by:
     * wi = -wo + 2 * dot(wo, wh) * wh
     */
    Color BlinnMaterial::sampleBSDF(const Fragment& fragment, 
        const Vector3& wo, const BSDFSample& bsdfSample, Vector3* wi, 
        float* pdf, BSDFType type, BSDFType* sampledType) const {
        BSDFType materialType = BSDFType(BSDFGlossy | BSDFReflection);
        if(!matchType(type, materialType)) {
            *pdf = 0.0f;
            return Color::Black;
        }
        float u1 = bsdfSample.uDirection[0];
        float u2 = bsdfSample.uDirection[1];
        float exp = mExp->lookup(fragment);
        float cosTheta = pow(u1, 1.0f / (exp + 1.0f));
        float sinTheta = sqrtf(max(0.0f, 1.0f - cosTheta * cosTheta));
        float phi = u2 * TWO_PI; 

        Vector3 whLocal(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
        // flip it if wo and normal at different sides of hemisphere
        if(dot(wo, fragment.getNormal()) < 0.0f) {
            whLocal *= -1.0f;
        }
        Matrix3 shadeToWorld = fragment.getWorldToShade().transpose();
        Vector3 wh = shadeToWorld * whLocal;;
        *wi = -wo + 2.0f * dot(wo, wh) * wh;
        *pdf = this->pdf(fragment, wo, *wi, materialType);
        if(sampledType) {
            *sampledType = materialType;
        }
        return bsdf(fragment, wo, *wi, materialType);
    }

    /* 
     * since dwh = dwo / (4 * dot(wo, wh))
     * pdf(wo) = pdf(wh) / (4 * dot(wo, wh))
     */
    float BlinnMaterial::pdf(const Fragment& fragment,
        const Vector3& wo, const Vector3& wi,
        BSDFType type) const {
        if(!matchType(type, BSDFType(BSDFGlossy | BSDFReflection))) {
            return 0.0f;
        }

        if(!sameHemisphere(fragment, wo, wi)) {
            return 0.0f;
        }
        Vector3 wh = normalize(wo + wi);
        float cosThetah = absdot(wh, fragment.getNormal());
        float exp = mExp->lookup(fragment);
        return (exp + 1.0f) * pow(cosThetah, exp) / 
            (TWO_PI * 4.0f * dot(wo, wh));
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
                    specularReflectDieletric(fragment, wo, wi, mEtai, mEtat);
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
            float reflect = specularReflectDieletric(fragment, wo, 
                &wReflect, mEtai, mEtat);
            float refract = specularRefract(fragment, wo, 
                &wRefract, mEtai, mEtat);
            // use fresnel factor to do importance sampling
            float fresnel = reflect * absdot(wReflect, fragment.getNormal());
            float reflectChance = fresnel;
            bool doReflect = mRNG->randomFloat() < reflectChance;

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


    Color MirrorMaterial::sampleBSDF(const Fragment& fragment, 
        const Vector3& wo, const BSDFSample& bsdfSample, Vector3* wi, 
        float* pdf, BSDFType type, BSDFType* sampledType) const {
        BSDFType materialType = BSDFType(BSDFSpecular | BSDFReflection);
        if(!matchType(type, materialType)) {
            *pdf = 0.0f;
            return Color::Black;
        }

        Color f = mReflectFactor->lookup(fragment) * 
            specularReflectConductor(fragment, wo, wi, mEta, mK);
        *pdf = 1.0f;
        if(sampledType) {
            *sampledType = materialType;
        }
        return f;
    }
}
