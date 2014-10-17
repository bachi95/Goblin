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
        VolumeRegion(Color absorption, Color emission, Color scatter,
            float g, float stepSize, const BBox& b, const Transform& toWorld);

        Color getAbsorption(const Vector3& p) const;
        Color getEmission(const Vector3& p) const;
        Color getScatter(const Vector3& p) const;
        Color getAttenuation(const Vector3& p) const;
        float getSampleStepSize() const;
        bool intersect(const Ray& ray, float* tMin, float* tMax) const;
        Color opticalThickness(const Ray& ray) const;
        float phase(const Vector3& p, const Vector3& wi, 
            const Vector3& wo) const;
        Color transmittance(const Ray& ray) const;
    private:
        Color mAbsorption;
        Color mEmission;
        Color mScatter;
        // param for Henyey-Greenstein evaluation
        float mG;
        // sample step size for source term and optical thickness
        float mStepSize; 
        BBox mLocalRegion;
        Transform mToWorld;
    };

    inline VolumeRegion::VolumeRegion(Color absorption, Color emission, 
        Color scatter, float g, float stepSize, 
        const BBox& b, const Transform& toWorld): mAbsorption(absorption),
        mEmission(emission), mScatter(scatter), mG(g), mStepSize(stepSize),
        mLocalRegion(b), mToWorld(toWorld) {}

    inline Color VolumeRegion::getAbsorption(const Vector3& p) const {
        bool inside = mLocalRegion.contain(mToWorld.invertPoint(p));
        return inside ? mAbsorption : Color(0.0f);
    }

    inline Color VolumeRegion::getEmission(const Vector3& p) const {
        bool inside = mLocalRegion.contain(mToWorld.invertPoint(p));
        return inside ? mEmission : Color(0.0f);
    }

    inline Color VolumeRegion::getScatter(const Vector3& p) const {
        bool inside = mLocalRegion.contain(mToWorld.invertPoint(p));
        return inside ? mScatter : Color(0.0f);
    }

    inline Color VolumeRegion::getAttenuation(const Vector3& p) const {
        bool inside = mLocalRegion.contain(mToWorld.invertPoint(p));
        return inside ? mAbsorption + mScatter : Color(0.0f);
    }

    inline float VolumeRegion::getSampleStepSize() const {
        return mStepSize;
    }


    class ParamSet;

    class VolumeCreator : public 
        Creator<VolumeRegion, const ParamSet&> {
    public:
        VolumeRegion* create(const ParamSet& params) const;
    };



}

#endif //GOBLIN_VOLUME_H

