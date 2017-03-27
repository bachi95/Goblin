#include "GoblinVolume.h"
#include "GoblinParamSet.h"
#include "GoblinRay.h"

namespace Goblin {

    bool VolumeRegion::intersect(const Ray& ray, 
        float* tMin, float* tMax) const {
        return mLocalRegion.intersect(mToWorld.invertRay(ray), tMin, tMax);
    }

    Color VolumeRegion::opticalThickness(const Ray& ray) const {
        float tMin, tMax;
        if(!mLocalRegion.intersect(mToWorld.invertRay(ray), &tMin, &tMax)) {
            return Color(0.0f);
        }
        return length(ray(tMax) - ray(tMin)) * (mAbsorption + mScatter);
    }

    float VolumeRegion::phase(const Vector3& p, const Vector3& wi, 
        const Vector3& wo) const {
        if(!mLocalRegion.contain(mToWorld.invertPoint(p))) {
            return 0.0f;
        }
        return phaseHG(wi, wo, mG);
    }

    Color VolumeRegion::transmittance(const Ray& ray) const {
        Color tau = opticalThickness(ray);
        return Color(exp(-tau.r), exp(-tau.g), exp(-tau.b));
    }


    VolumeRegion* VolumeCreator::create(const ParamSet& params) const {
        Color absorption = params.getColor("absorption");
        Color emission = params.getColor("emission");
        Color scatter = params.getColor("scatter");
        float g = params.getFloat("g", 0.0f);
        float stepSize = params.getFloat("step_size", 0.1f);
        int sampleNum = params.getInt("sample_num", 5);
        Vector3 vMin = params.getVector3("box_min");
        Vector3 vMax = params.getVector3("box_max");
        BBox b(vMin, vMax);

        Vector3 position = params.getVector3("position");
        Vector4 v = params.getVector4("orientation", 
            Vector4(1.0f, 0.0f, 0.0f, 0.0f));
        Quaternion orientation(v[0], v[1], v[2], v[3]);
        Transform toWorld(position, orientation, Vector3(1.0f, 1.0f, 1.0f));

        return new VolumeRegion(absorption, emission, scatter, g, stepSize,
            sampleNum, b, toWorld);
    }

}
