#include "GoblinMaterial.h"
#include "GoblinVector.h"
#include "GoblinMatrix.h"
#include "GoblinVertex.h"
#include "GoblinGeometry.h"
#include "GoblinLight.h"
#include "GoblinParamSet.h"
#include "GoblinRay.h"
#include "GoblinSampler.h"
#include "GoblinScene.h"
#include "GoblinVolume.h"

namespace Goblin {

    BSDFSampleIndex::BSDFSampleIndex(SampleQuota* sampleQuota, 
        int requestNum) {
        SampleIndex oneDIndex = sampleQuota->requestOneDQuota(requestNum);
        SampleIndex twoDIndex = sampleQuota->requestTwoDQuota(requestNum);
        // theoretically this two should be the same...
        // this is just a paranoid double check
        samplesNum = min(oneDIndex.sampleNum, twoDIndex.sampleNum);
        componentIndex = oneDIndex.offset;
        directionIndex = twoDIndex.offset;
    }

    BSDFSample::BSDFSample(const RNG& rng) {
        uComponent = rng.randomFloat();
        uDirection[0] = rng.randomFloat();
        uDirection[1] = rng.randomFloat();
    }

    BSDFSample::BSDFSample(const Sample& sample,
        const BSDFSampleIndex& index, uint32_t n) {
        uComponent = sample.u1D[index.componentIndex][n];
        uDirection[0] = sample.u2D[index.directionIndex][2 * n];
        uDirection[1] = sample.u2D[index.directionIndex][2 * n + 1];
    }

    BSSRDF::BSSRDF(const ColorTexturePtr& absorb, 
        const ColorTexturePtr& scatterPrime, float eta, float g):
        mAbsorb(absorb), mScatterPrime(scatterPrime), 
        mEta(eta), mG(g) {
        float fdr = Fdr(eta);
        mA = (1.0f + fdr) / (1.0f - fdr);
    }


    BSSRDF::BSSRDF(const Color& Kd, const Color& diffuseMeanFreePath, 
        float eta, float g): mEta(eta), mG(g) {
        float fdr = Fdr(eta);
        mA = (1.0f + fdr) / (1.0f - fdr);
        Color absorb, scatterPrime;
        convertFromDiffuse(Kd, diffuseMeanFreePath, mA, 
            &absorb, &scatterPrime); 
        mAbsorb = ColorTexturePtr(new ConstantTexture<Color>(absorb));
        mScatterPrime = 
            ColorTexturePtr(new ConstantTexture<Color>(scatterPrime));
    }

    Color BSSRDF::Rd(const Fragment& fragment, float d2) const {
        // see Donner. C 2006 Chapter 5 for the full derivation 
        // of the following disffusion dipole approximation equation
        Color sigmaA = mAbsorb->lookup(fragment);
        Color sigmaSPrime = mScatterPrime->lookup(fragment);
        Color sigmaTPrime = sigmaA + sigmaSPrime;
        Color sigmaTr = sqrtColor(3.0f * sigmaA * sigmaTPrime);
        Color one(1.0f);
        Color zr = one / sigmaTPrime;
        // zv = zr + 4AD where D = 1/(3 * sigmaT') = zr / 3
        Color zv = zr * (1.0f + 4.0f / 3.0f * mA);
        Color dr = sqrtColor(zr * zr + Color(d2));
        Color dv = sqrtColor(zv * zv + Color(d2));

        Color alphaPrime = sigmaSPrime / sigmaTPrime;
        Color sTrDr = sigmaTr * dr;
        Color sTrDv = sigmaTr * dv;
        Color rd = 0.25f * INV_PI * alphaPrime * (
            (zr * (one + sTrDr) * expColor(-sTrDr) / (dr * dr * dr)) +
            (zv * (one + sTrDv) * expColor(-sTrDv) / (dv * dv * dv)));
        return clampColor(rd);
    }

    float BSSRDF::MISWeight(const Fragment& fo, const Fragment& fi,
        BSSRDFSampleAxis mainAxis, float pdf, float sigmaTr,
        float Rmax) const {
        float weight = 0.0f;
        // The U, V, N ratio is 1 : 1 : 2, and we do MIS with
        // power heuristic, so 1: 1: 4 
        const Vector3& pwo = fo.getPosition();
        const Vector3& pwi = fi.getPosition();
        const Vector3& ni = fi.getNormal();
        if (mainAxis == NAxis) {
            const Vector3& u = normalize(fo.getDPDU());
            const Vector3& v = normalize(fo.getDPDV());
            float UPdf = 0.25f * 
                gaussianSample2DPdf(pwo, pwi, u, sigmaTr, Rmax) * 
                absdot(u, ni);
            float VPdf = 0.25f * 
                gaussianSample2DPdf(pwo, pwi, v, sigmaTr, Rmax) *
                absdot(v, ni);
            float numerator = 4 *  pdf * pdf;
            weight = numerator / (numerator + UPdf* UPdf + VPdf * VPdf);
        } else if (mainAxis == UAxis) {
            const Vector3& n = fo.getNormal();
            const Vector3& v = normalize(fo.getDPDV());
            float NPdf = 0.5f * 
                gaussianSample2DPdf(pwo, pwi, n, sigmaTr, Rmax) * 
                absdot(n, ni);
            float VPdf = 0.25f * 
                gaussianSample2DPdf(pwo, pwi, v, sigmaTr, Rmax) *
                absdot(v, ni);
            float numerator = pdf * pdf;
            weight = numerator / (4 * NPdf * NPdf + numerator + VPdf * VPdf);
        } else if (mainAxis == VAxis) {
            const Vector3& n = fo.getNormal();
            const Vector3& u = normalize(fo.getDPDU());
            float NPdf = 0.5f * 
                gaussianSample2DPdf(pwo, pwi, n, sigmaTr, Rmax) * 
                absdot(n, ni);
            float UPdf = 0.25f * 
                gaussianSample2DPdf(pwo, pwi, u, sigmaTr, Rmax) *
                absdot(u, ni);
            float numerator = pdf * pdf;
            weight = numerator / (4 * NPdf * NPdf + UPdf * UPdf + numerator);
        }
        return weight;
    }

    BSSRDFSampleAxis BSSRDF::sampleProbeRay(const Fragment& fragment,
        const BSSRDFSample& sample, float sigmaTr, float Rmax,
        Ray* probeRay, float* pdf) const {
        Matrix3 shadeToWorld = fragment.getWorldToShade().transpose();
        const Vector3& pwo = fragment.getPosition();
        // sample disk with gaussian falloff
        Vector2 pSample = gaussianSample2D(
            sample.uDisc[0], sample.uDisc[1], sigmaTr, Rmax);
        float halfProbeLength = sqrt(Rmax * Rmax - pSample.squaredLength());
        // figure out we should sample alone U or V or N 
        // The chance to pick up U, V, N is 1 : 1 : 2
        BSSRDFSampleAxis axis;
        if (sample.uPickAxis <= 0.5f) {
            probeRay->o =  pwo + shadeToWorld * 
                Vector3(pSample.x, pSample.y, -halfProbeLength);
            probeRay->d = fragment.getNormal();
            axis = NAxis;
            *pdf = 0.5f;
        } else if (sample.uPickAxis <= 0.75f) {
            probeRay->o =  pwo + shadeToWorld * 
                Vector3(-halfProbeLength, pSample.x, pSample.y);
            probeRay->d = normalize(fragment.getDPDU());
            axis = UAxis;
            *pdf = 0.25f;
        } else {
            probeRay->o =  pwo + shadeToWorld * 
                Vector3(pSample.y, -halfProbeLength, pSample.x);
            probeRay->d = normalize(fragment.getDPDV());
            axis = VAxis;
            *pdf = 0.25f;
        }
        probeRay->mint = 0.0f;
        probeRay->maxt = 2.0f * halfProbeLength;
        *pdf *= gaussianSample2DPdf(pSample.x, pSample.y, sigmaTr, Rmax); 
        return axis;
    }

    Color BSSRDF::getSigmaTr(const Fragment& fragment) const {
        Color sigmaA = mAbsorb->lookup(fragment);
        Color sigmaSPrime = mScatterPrime->lookup(fragment);
        Color sigmaTPrime = sigmaA + sigmaSPrime;
        return sqrtColor(3.0f * sigmaA * sigmaTPrime);
    }

    float BSSRDF::phase(const Vector3& wi, const Vector3& wo) const {
        return phaseHG(wi, wo, mG);
    }

    void BSSRDF::convertFromDiffuse(const Color& Kd, 
        const Color& diffuseMeanFreePath, float A, 
        Color* absorb, Color* scatterPrime) {
        float diffuse[3] = {Kd.r, Kd.g, Kd.b};
        float sigmaTr[3] = {
            1.0f / diffuseMeanFreePath.r,
            1.0f / diffuseMeanFreePath.g, 
            1.0f / diffuseMeanFreePath.b};

        float sigmaSPrime[3];
        float sigmaA[3];
        // RGB channel
        for (int i = 0; i < 3; ++i) {
            // use few binary search iterations to approximate alpha prime
            // since there is no easy way to inverse diffuseReflectance method
            // directly...
            float alphaLow = 0.0f;
            float alphaHigh = 1.0f;
            for (int j = 0; j < 16; ++j) {
                float alphaMid =  0.5f * (alphaLow + alphaHigh);                
                float rdMid = diffuseReflectance(alphaMid, A);
                if (rdMid > diffuse[i]) {
                    alphaHigh = alphaMid;
                } else {
                    alphaLow = alphaMid;
                }
            }
            float alphaPrime = 0.5f * (alphaLow + alphaHigh); 
            float sigmaTPrime = sigmaTr[i] / 
                sqrt(3.0f * (1.0f - alphaPrime));
            sigmaSPrime[i] = alphaPrime * sigmaTPrime;
            sigmaA[i] = sigmaTPrime - sigmaSPrime[i];
        }
        *scatterPrime = Color(sigmaSPrime[0], sigmaSPrime[1], sigmaSPrime[2]);
        *absorb = Color(sigmaA[0], sigmaA[1], sigmaA[2]);
    }

    float BSSRDF::diffuseReflectance(float alphaPrime, float A) {
        float sqrtTerm = sqrt(3.0f * (1.0f - alphaPrime));
        return 0.5f * alphaPrime * 
            (1.0f + exp(-(4.0f / 3.0f) * A * sqrtTerm)) *
            exp(-sqrtTerm);
    }

    void BumpShaders::evaluate(Fragment* fragment) const {
        /*
         * let bumpmap(u, v) be D
         * p'(u, v) = p(u, v) + n(u, v) * D
         * we want to get the new normal n'(u, v),
         * which is cross(dp'du, dp'dv)
         * dp'du = d(p'(u, v)/du = dp/du + dD/du * n + dn/du * D
         * we have dp/du already
         * dn/du * D can be ignored when D is small
         * dD/du can be approximated by forward differential approximation
         * (D(u + du, v) - D(u, v)) / du
         * dp'dv is calculated with the same way above
         */
        if (bumpMap) {
            Vector3 p = fragment->getPosition();
            Vector3 n = fragment->getNormal();
            Vector2 uv = fragment->getUV();
            float bumpD = bumpMap->lookup(*fragment);

            float du = 0.002f;
            Fragment fdu = *fragment; 
            fdu.setPosition(p + du * fragment->getDPDU());
            fdu.setUV(uv + Vector2(du, 0.0f));
            float bumpDdu = bumpMap->lookup(fdu);
            Vector3 bumpDPDU = fragment->getDPDU() + 
                (bumpDdu - bumpD) / du * n;

            float dv = 0.002f;
            Fragment fdv = *fragment;
            fdv.setPosition(p + dv * fragment->getDPDV());
            fdv.setUV(uv + Vector2(0.0f, dv));
            float bumpDdv = bumpMap->lookup(fdv);
            Vector3 bumpDPDV = fragment->getDPDV() +
                (bumpDdv - bumpD) / dv * n;

            Vector3 bumpN = normalize(cross(bumpDPDU, bumpDPDV));
            // case that intersect back face
            if (dot(bumpN, n) < 0.0f) {
                bumpN *= -1.0f;
            }
            // now set back the bumped value
            fragment->setNormal(bumpN);
            fragment->setDPDU(bumpDPDU);
            fragment->setDPDV(bumpDPDV);
        }

        if (normalMap) {
            // Decode normal map
            Color colorN = normalMap->lookup(*fragment);
            Vector3 nShade(colorN.r, colorN.g, colorN.b); 
            nShade = 2.0f * nShade - Vector3(1.0f, 1.0f, 1.0f);
            // mathmatically it should be shadeToWorld.inverse().transpose()
            // but since getWorldToShade guarantee to deliver orthognal matrix
            // it's fine that we directly use its transpose to transform normal
            Matrix3 shadeToWorld = fragment->getWorldToShade().transpose();
            Vector3 nWorld = normalize(shadeToWorld * nShade);
            if (dot(nWorld, fragment->getNormal()) < 0.0f) {
                nWorld *= -1.0f;
            }
            fragment->setNormal(nWorld);
        }
        return;
    }

    BSDFType Material::getSampleType(const Vector3& wo, const Vector3& wi,
        const Fragment& fragment, BSDFType type) const {
        const Vector3& n = fragment.getNormal();
        if (dot(n, wo) * dot(n, wi) > 0.0f) {
            type = BSDFType(type & ~BSDFTransmission);
        } else {
            type = BSDFType(type & ~BSDFReflection);
        }
        return type;
    }

    void Material::perturb(Fragment* fragment) const {
        mBumpShaders.evaluate(fragment);
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
        if (!entering) {
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
        if (cosi <= 0.0f) {
            return 0.0f;
        }
        float f = fresnelConductor(cosi, eta, k);
        *wi = 2 * cosi * n - wo;
        float cosr = cosi;
        return f / cosr;
    }

    float Material::specularRefract(const Fragment& fragment,
        const Vector3& wo, Vector3* wi, float etao, float etai,
        BSDFMode mode) const {
        // wo(face outward) is the input ray(with n form angle i), 
        // wi is the refract ray(with -n form angle t)
        Vector3 n = fragment.getNormal();
        // wo(face outward) is the input ray(with n form angle i), 
        // wi is the refract ray
        float coso = dot(n, wo);
        float et = etao;
        float ei = etai;
        bool entering = coso > 0.0f;
        if (!entering) {
            swap(ei, et);
            n = -n;
            coso = -coso;
        }
        float f = fresnelDieletric(coso, et, ei);
        // total reflection
        if (f == 1.0f) {
            return 0.0f;
        }
        /*
         * Wi = -N * cosi - WoPerpN * sini / sino =
         * -N * cosi - (sini / sino) * (Wo - dot(N, Wo) * N) =
         * -N * sqrt(1 - (etao / etai)^2 * (1 - (dot(N, Wo))^2))) +
         * etao / etai * (Wo - dot(N, Wo) * N) =
         * N * (etao / etai * dot(N, Wo) -
         * sqrt(1 - (etao / etai)^2(1 - (dot(N, Wo))^2)) -
         * etao / etai * Wo
         */
        float eta = et / ei;
        *wi = normalize(n * (eta * coso -
            sqrt(max(0.0f, 1.0f - eta * eta * (1.0f - coso * coso)))) -
            eta * wo);
        /*
         * see Veach 97 Chapter 5: The Sources of Non-Symmetric Scattering
         * for detail derivation, solid angle would be "squeezed" from
         * smaller IOR to bigger IOR and cause radiance increase, while
         * it doesn't apply to importance
         */
        return mode == BSDFRadiance ?
            eta * eta * (1.0f - f) / absdot(*wi, n) :
            (1.0f - f) / absdot(*wi, n);
    }

    // close approximation fresnel for dieletric material
    float Material::fresnelDieletric(float cosi, float etai, 
        float etat) {
        cosi = clamp(cosi, -1.0f, 1.0f);
        float sint = (etai / etat) * sqrt(max(0.0f, 1.0f - cosi * cosi));
        // total reflection
        if (sint >= 1.0f) {
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

    float Material::fresnelConductor(float cosi, float eta, float k) {
        float tmp = (eta * eta + k * k);
        float cosi2 = cosi * cosi;
        float rParl2 = (tmp * cosi2 - 2.0f * eta * cosi + 1.0f) /
            (tmp * cosi2 + 2.0f * eta * cosi + 1.0f);
        float rPerp2 = (tmp - 2.0f * eta * cosi + cosi2) /
            (tmp + 2.0f * eta * cosi + cosi2);
        return (rParl2 + rPerp2) * 0.5f;
    }

    Vector3 specularRefract(const Vector3& wo, const Vector3& n,
        float etai, float etat) {
        /*
         * Wi = -N * cost - WoPerpN * sint / sini =
         * -N * cost - (sint / sini) * (Wo - dot(N, Wo) * N) =
         * -N * sqrt(1 - (etai / etat)^2 * (1 - (dot(N, Wo))^2))) +
         * etai / etat * (Wo - dot(N, Wo) * N) =
         * N * (etai/etat * dot(N, Wo) - 
         * sqrt(1 - (etai/etat)^2(1 - (dot(N, Wo))^2)) -
         * etai / etat * Wo
         */
        float eta = etai / etat;
        float cosi = absdot(n, wo);
        return normalize(n * (eta * cosi - 
            sqrt(max(0.0f, 1.0f - eta * eta * (1.0f - cosi * cosi)))) - 
            eta * wo);
    }


    Color LambertMaterial::bsdf(const Fragment& fragment, const Vector3& wo, 
        const Vector3& wi, BSDFType type, BSDFMode mode) const {
        Color f(Color::Black);
        type = getSampleType(wo, wi, fragment, type);
        if (matchType(type, getType())) {
            f += mDiffuseFactor->lookup(fragment) * INV_PI;
        }
        return f;
    }

    Color LambertMaterial::sampleBSDF(const Fragment& fragment, 
        const Vector3& wo, const BSDFSample& bsdfSample, Vector3* wi, 
        float* pdf, BSDFType type, BSDFType* sampledType,
        BSDFMode mode) const {
        BSDFType materialType = getType();
        if (!matchType(type, materialType)) {
            *pdf = 0.0f;
            return Color::Black;
        }
        float u1 = bsdfSample.uDirection[0];
        float u2 = bsdfSample.uDirection[1];
        Vector3 wiLocal = cosineSampleHemisphere(u1, u2);
        // flip it if wo and normal at different sides of hemisphere
        if (dot(wo, fragment.getNormal()) < 0.0f) {
            wiLocal *= -1.0f;
        }
        Matrix3 shadeToWorld = fragment.getWorldToShade().transpose();
        *wi = shadeToWorld * wiLocal;;
        *pdf = this->pdf(fragment, wo, *wi, materialType);
        if (sampledType) {
            *sampledType = materialType;
        }
        return mDiffuseFactor->lookup(fragment) * INV_PI;
    }

    float LambertMaterial::pdf(const Fragment& fragment,
        const Vector3& wo, const Vector3& wi,
        BSDFType type) const {
        if (!matchType(type, getType())) {
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
        const Vector3& wi, BSDFType type, BSDFMode mode) const {
        type = getSampleType(wo, wi, fragment, type);
        if (matchType(type, getType())) {
            Vector3 n = fragment.getNormal();
            float cosi = absdot(n, wi);
            float coso = absdot(n, wo);
            if (cosi == 0.0f || coso == 0.0f) {
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
            if (mFresnelType == Dieletric) {
                F = fresnelDieletric(woDotWh, 1.0f, mEta);
            } else if (mFresnelType == Conductor) {
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
        float* pdf, BSDFType type, BSDFType* sampledType,
        BSDFMode mode) const {
        BSDFType materialType = getType();
        if (!matchType(type, materialType)) {
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
        if (dot(wo, fragment.getNormal()) < 0.0f) {
            whLocal *= -1.0f;
        }
        Matrix3 shadeToWorld = fragment.getWorldToShade().transpose();
        Vector3 wh = shadeToWorld * whLocal;;
        *wi = -wo + 2.0f * dot(wo, wh) * wh;
        *pdf = this->pdf(fragment, wo, *wi, materialType);
        if (sampledType) {
            *sampledType = materialType;
        }
        return bsdf(fragment, wo, *wi, materialType, mode);
    }

    /* 
     * since dwh = dwo / (4 * dot(wo, wh))
     * pdf(wo) = pdf(wh) / (4 * dot(wo, wh))
     */
    float BlinnMaterial::pdf(const Fragment& fragment,
        const Vector3& wo, const Vector3& wi,
        BSDFType type) const {
        if (!matchType(type, getType())) {
            return 0.0f;
        }

        if (!sameHemisphere(fragment, wo, wi)) {
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
        float* pdf, BSDFType type, BSDFType* sampledType,
        BSDFMode mode) const {
        int nMatch = 0;
        if (matchType(type, BSDFType(BSDFSpecular | BSDFReflection))) {
            nMatch++;
        }
        if (matchType(type, BSDFType(BSDFSpecular | BSDFTransmission))) {
            nMatch++;
        }

        Color f(Color::Black);
        if (nMatch == 1) {
            if (matchType(type, BSDFReflection)) {
                f = mReflectFactor->lookup(fragment) * 
                    specularReflectDieletric(fragment, wo, wi, mEtai, mEtat);
                if (sampledType) {
                    *sampledType = BSDFType(BSDFSpecular | BSDFReflection);
                }
            } else {
                f = mRefractFactor->lookup(fragment) *
                    specularRefract(fragment, wo, wi, mEtai, mEtat, mode);
                if (sampledType) {
                    *sampledType = BSDFType(BSDFSpecular | BSDFTransmission);
                }
            }
            *pdf = 1.0f;
        } else if (nMatch == 2) {
            Vector3 wReflect;
            Vector3 wRefract;
            float reflect = specularReflectDieletric(fragment, wo, 
                &wReflect, mEtai, mEtat);
            float refract = specularRefract(fragment, wo, 
                &wRefract, mEtai, mEtat, mode);
            // use fresnel factor to do importance sampling
            float fresnel = reflect * absdot(wReflect, fragment.getNormal());
            float reflectChance = fresnel;
            bool doReflect = bsdfSample.uComponent < reflectChance;

            if (doReflect) {
                f = mReflectFactor->lookup(fragment) * reflect;
                *wi = wReflect;
                if (sampledType) {
                    *sampledType = BSDFType(BSDFSpecular | BSDFReflection);
                }
                *pdf = reflectChance;
            } else {
                f = mRefractFactor->lookup(fragment) * refract;
                *wi = wRefract;
                if (sampledType) {
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
        float* pdf, BSDFType type, BSDFType* sampledType,
        BSDFMode mode) const {
        BSDFType materialType = getType();
        if (!matchType(type, materialType)) {
            *pdf = 0.0f;
            return Color::Black;
        }

        Color f = mReflectFactor->lookup(fragment) * 
            specularReflectConductor(fragment, wo, wi, mEta, mK);
        *pdf = 1.0f;
        if (sampledType) {
            *sampledType = materialType;
        }
        return f;
    }

    Color SubsurfaceMaterial::sampleBSDF(const Fragment& fragment, 
        const Vector3& wo, const BSDFSample& bsdfSample, Vector3* wi, 
        float* pdf, BSDFType type, BSDFType* sampledType,
        BSDFMode mode) const {
        BSDFType materialType = getType();
        if (!matchType(type, materialType)) {
            *pdf = 0.0f;
            return Color::Black;
        }

        Color f = mReflectFactor->lookup(fragment) * 
            specularReflectDieletric(fragment, wo, wi, 1.0f, mEta);
        *pdf = 1.0f;
        if (sampledType) {
            *sampledType = materialType;
        }
        return f;
    }

    Color MaskMaterial::bsdf(const Fragment& fragment, const Vector3& wo,
        const Vector3& wi, BSDFType type, BSDFMode mode) const {
        if (type == BSDFNull) {
            return Color::Black;
        } else {
            float alpha = mAlphaMask->lookup(fragment);
            return alpha * mMaskedMaterial->bsdf(fragment, wo, wi, type);
        }
    }

    Color MaskMaterial::sampleBSDF(const Fragment& fragment, 
        const Vector3& wo, const BSDFSample& bsdfSample, Vector3* wi,
        float* pdf, BSDFType type, BSDFType* sampledType,
        BSDFMode mode) const {
        bool sampleAlpha = matchType(type, BSDFNull);
        // if this request ask us to evaluate masked material
        bool sampleMasked = type != BSDFNull;

        Color result(0.0f);
        *pdf = 0.0f;
        float alpha = mAlphaMask->lookup(fragment);
        if (sampleMasked && sampleAlpha) {
            float maskedProb = alpha;
            if (bsdfSample.uComponent < maskedProb) {
                result = alpha * mMaskedMaterial->sampleBSDF(fragment, wo, 
                    bsdfSample, wi, pdf, type, sampledType);
                *pdf *= maskedProb;
            } else {
                result = (1.0f - alpha) * mTransparentColor->lookup(fragment);
                *wi = -normalize(wo);
                *pdf = 1.0f - maskedProb;
                if (sampledType) {
                    *sampledType = BSDFNull;
                }
            }
        } else if (sampleMasked) {
            result = alpha * mMaskedMaterial->sampleBSDF(fragment, wo, 
                bsdfSample, wi, pdf, type, sampledType);
        } else if (sampleAlpha) {
            result = (1.0f - alpha) * mTransparentColor->lookup(fragment);
            *wi = -normalize(wo);
            *pdf = 1.0f;
            if (sampledType) {
                *sampledType = BSDFNull;
            }
        }
        return result;
    }

    float MaskMaterial::pdf(const Fragment& fragment, const Vector3& wo, 
        const Vector3& wi, BSDFType type) const {
        bool sampleAlpha = matchType(type, BSDFNull);
        // if this request ask us to evaluate masked material
        bool sampleMasked = type != BSDFNull;
        float pdf = 0.0f;
        if (sampleMasked && sampleAlpha) {
            pdf = mAlphaMask->lookup(fragment) * 
                mMaskedMaterial->pdf(fragment, wo, wi, type);
        } else if (sampleMasked) {
            pdf = mMaskedMaterial->pdf(fragment, wo, wi, type);
        } else if (sampleAlpha) {
            pdf = 0.0f;
        }
        return pdf;
    }

    static BumpShaders getBumpShaders(const ParamSet& params,
        const SceneCache& sceneCache) {
        FloatTexturePtr bump;
        if (params.hasString("bumpmap")) {
            bump = sceneCache.getFloatTexture(params.getString("bumpmap"));
        }
        ColorTexturePtr normal;
        if (params.hasString("normalmap")) {
            normal = sceneCache.getColorTexture(params.getString("normalmap"));
        }
        return BumpShaders(bump, normal);
    }

    Material* LambertMaterialCreator::create(const ParamSet& params,
        const SceneCache& sceneCache) const {
        string textureName = params.getString("Kd");
        ColorTexturePtr Kd = sceneCache.getColorTexture(textureName);
        BumpShaders bump = getBumpShaders(params, sceneCache);
        return new LambertMaterial(Kd, bump);
    }


    Material* BlinnMaterialCreator::create(const ParamSet& params,
        const SceneCache& sceneCache) const {
        string glossyTextureName = params.getString("Kg");
        string expTextureName = params.getString("exponent");
        ColorTexturePtr Kg = 
            sceneCache.getColorTexture(glossyTextureName);
        FloatTexturePtr exp = 
            sceneCache.getFloatTexture(expTextureName);
        float index = params.getFloat("index", 1.5f);
        float absorption = params.getFloat("k", -1.0f);
        BumpShaders bump = getBumpShaders(params, sceneCache);
        Material* material;
        if (absorption > 0.0f) {
            // conductor
            material = new BlinnMaterial(Kg, exp, index, absorption, bump);
        } else {
            // dieletric
            material = new BlinnMaterial(Kg, exp, index, bump);
        }
        return material;
    }


    Material* TransparentMaterialCreator::create(const ParamSet& params,
        const SceneCache& sceneCache) const {
        string reflectTextureName = params.getString("Kr");
        string refractTextureName = params.getString("Kt");
        ColorTexturePtr Kr = 
            sceneCache.getColorTexture(reflectTextureName);
        ColorTexturePtr Kt = 
            sceneCache.getColorTexture(refractTextureName);
        float index = params.getFloat("index", 1.5f);
        BumpShaders bump = getBumpShaders(params, sceneCache);
        return new TransparentMaterial(Kr, Kt, index, bump);
    }


    Material* MirrorMaterialCreator::create(const ParamSet& params,
        const SceneCache& sceneCache) const {
        string reflectTextureName = params.getString("Kr");
        ColorTexturePtr Kr = 
            sceneCache.getColorTexture(reflectTextureName);
        // use aluminium by default
        float index = params.getFloat("index", 0.8f);
        float absorption = params.getFloat("k", 6.0f);
        BumpShaders bump = getBumpShaders(params, sceneCache);
        return new MirrorMaterial(Kr, index, absorption, bump);
    }


    Material* SubsurfaceMaterialCreator::create(const ParamSet& params,
        const SceneCache& sceneCache) const {

        float index = params.getFloat("index", 1.5f);
        float g = params.getFloat("g", 0.0f);
        ColorTexturePtr Kr;
        if (params.hasString("Kr")) {
            Kr = sceneCache.getColorTexture(params.getString("Kr"));
        } else {
            Kr = ColorTexturePtr(
                new ConstantTexture<Color>(Color(1.0f)));
        }
        BumpShaders bump = getBumpShaders(params, sceneCache);

        if (params.hasColor("Kd")) {
            Color Kd = params.getColor("Kd");
            Color diffuseMeanFreePath = params.getColor("mean_free_path", 
                Color(1.0f));
            return new SubsurfaceMaterial(Kd, diffuseMeanFreePath, Kr, 
                index, g, bump);

        } else {
            ColorTexturePtr absorb;
            if (params.hasString("absorb")) {
                absorb = sceneCache.getColorTexture(
                    params.getString("absorb"));
            } else {
                Color marbleAbsorb(0.0021f, 0.0041f, 0.0071f);
                absorb = ColorTexturePtr(
                    new ConstantTexture<Color>(marbleAbsorb));
            }

            ColorTexturePtr scatterPrime;
            if (params.hasString("scatter_prime")) {
                scatterPrime = sceneCache.getColorTexture(
                    params.getString("scatter_prime"));
            } else {
                Color marbleScatterp(2.19f, 2.62f, 3.00f);
                scatterPrime = ColorTexturePtr(
                    new ConstantTexture<Color>(marbleScatterp));
            }

            return new SubsurfaceMaterial(absorb, scatterPrime, Kr, 
                index, g, bump);
        }
    }


    Material* MaskMaterialCreator::create(const ParamSet& params,
        const SceneCache& sceneCache) const {
        FloatTexturePtr alphaMask;
        if (params.hasString("alpha")) {
            alphaMask = sceneCache.getFloatTexture(params.getString("alpha"));
        } else {
            alphaMask = FloatTexturePtr(new ConstantTexture<float>(1.0f));
            cout << "no feed in alpha" << endl;
        }
        ColorTexturePtr transparentColor;
        if (params.hasString("transparent_color")) {
            transparentColor= sceneCache.getColorTexture(
                params.getString("transparent_color"));
        } else {
            transparentColor = ColorTexturePtr(
                new ConstantTexture<Color>(Color::White));
            cout << "no feed in tr color" << endl;
        }
        MaterialPtr material = 
            sceneCache.getMaterial(params.getString("material"));
        return new MaskMaterial(alphaMask, transparentColor, material);
    }
}
