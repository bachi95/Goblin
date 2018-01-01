#ifndef GOBLIN_VOLUME_H
#define GOBLIN_VOLUME_H

#include "GoblinBBox.h"
#include "GoblinColor.h"
#include "GoblinFactory.h"
#include "GoblinTransform.h"

namespace Goblin {
    class Ray;

    class VolumeRegion {
    public:
        VolumeRegion(float g, float stepSize, int sampleNum, const BBox& b,
            const Transform& toWorld) :
            mG(g), mStepSize(stepSize), mSampleNum(sampleNum),
            mLocalRegion(b), mToWorld(toWorld) {}

        virtual ~VolumeRegion() {}

        // query absorbtion coefficient(sigmaA) at world space position p
        virtual Color getAbsorption(const Vector3& p) const = 0;

        // query emission coefficient(Lve) at world space position p
        virtual Color getEmission(const Vector3& p) const = 0;

        // query scatter coefficient(sigmaS) at world space position p
        virtual Color getScatter(const Vector3& p) const = 0;

        // query extinction coefficient(sigmaT) at world space position p
        // sigmaT = sigmaA + sigmaS
        virtual Color getAttenuation(const Vector3& p) const = 0;

        // calculate opticalThickness alone input ray
        // (the integration of extinction coefficient alone the ray)
        virtual Color opticalThickness(const Ray& ray) const = 0;

        // query ray marching step size (measured in world space)
        float getSampleStepSize() const {
            return mStepSize;
        }

        // query single scattering light integration sample count alone a ray
        int getLightSampleNum() const {
            return mSampleNum;
        }

        bool intersect(const Ray& ray, float* tMin, float* tMax) const;

        // evaluate phase function
        float phase(const Vector3& p, const Vector3& wi, 
            const Vector3& wo) const;

        // calculate transmittance alone input ray
        Color transmittance(const Ray& ray) const;

    protected:
        // param for Henyey-Greenstein evaluation
        float mG;
        // sample step size for source term and optical thickness
        float mStepSize;
        int mSampleNum;
        BBox mLocalRegion;
        Transform mToWorld;
    };

    // homogeneous volume implementation that absorb/emission/scatter
    // coefficients are constant
    class HomogeneousVolumeRegion : public VolumeRegion {
    public:
        HomogeneousVolumeRegion(const Color& absorption,
            const Color& emission, const Color& scatter,
            float g, float stepSize, int sampleNum, const BBox& b,
            const Transform& toWorld): VolumeRegion(g, stepSize, sampleNum, b, toWorld),
            mAbsorption(absorption), mEmission(emission), mScatter(scatter) {}

        Color getAbsorption(const Vector3& p) const {
            bool inside = mLocalRegion.contain(mToWorld.invertPoint(p));
            return inside ? mAbsorption : Color(0.0f);
        }

        Color getEmission(const Vector3& p) const {
            bool inside = mLocalRegion.contain(mToWorld.invertPoint(p));
            return inside ? mEmission : Color(0.0f);
        }

        Color getScatter(const Vector3& p) const {
            bool inside = mLocalRegion.contain(mToWorld.invertPoint(p));
            return inside ? mScatter : Color(0.0f);
        }

        Color getAttenuation(const Vector3& p) const {
            bool inside = mLocalRegion.contain(mToWorld.invertPoint(p));
            return inside ? mAbsorption + mScatter : Color(0.0f);
        }

        Color opticalThickness(const Ray& ray) const;

    private:
        Color mAbsorption;
        Color mEmission;
        Color mScatter;
    };

    // Henyey-Greenstein phase function.
    inline float phaseHG(const Vector3& wi, const Vector3& wo, float g) {
        if(g < 1e-3) {
            return 0.25f * INV_PI;
        } else {
            float cosTheta = dot(wi, wo);
            return 0.25f * INV_PI * (1.0f - g * g) / 
                powf(1.0f + g * g - 2.0f * g * cosTheta, 1.5f);
        }
    }


    class ParamSet;

    class HomogeneousVolumeCreator : public
        Creator<VolumeRegion, const ParamSet&> {
    public:
        VolumeRegion* create(const ParamSet& params) const;
    };

}

#endif //GOBLIN_VOLUME_H

