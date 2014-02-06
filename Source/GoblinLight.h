#ifndef GOBLIN_LIGHT_H
#define GOBLIN_LIGHT_H

#include "GoblinColor.h"
#include "GoblinVector.h"
#include "GoblinParamSet.h"
#include "GoblinUtils.h"
#include "GoblinGeometry.h"
#include "GoblinTransform.h"

#include <vector>

namespace Goblin {
    class Ray;
    class Quaternion;
    class CDF1D;

    class Light {
    public:
        enum Type {
            Point = 0,
            Directional = 1,
            Area = 2
        };
        virtual ~Light() {};
        virtual Color sampleL(const Vector3& p, float epsilon, 
            Vector3* wi, Ray* shadowRay) const = 0;
        ParamSet& getParams();
    protected:
        ParamSet mParams;
    };

    inline ParamSet& Light::getParams() {
        return mParams;
    }

    class PointLight : public Light {
    public:
        PointLight(const Color& intensity, const Vector3& position);
        Color sampleL(const Vector3& p, float epsilon, Vector3* wi,
            Ray* shadowRay) const;
        Color intensity;
        Vector3 position;
    };


    class DirectionalLight : public Light {
    public:
        DirectionalLight(const Color& radiance, const Vector3& direction);
        Color sampleL(const Vector3& p, float epsilon, Vector3* wi,
            Ray* shadowRay) const;
        Color radiance;
        Vector3 direction;
    };


    class GeometrySet {
    public:
        GeometrySet(const GeometryPtr& geometry);
        ~GeometrySet();
    private:
        GeometryList mGeometries;
        vector<float> mGeometriesArea;
        float mSumArea;
        CDF1D* mAreaDistribution;
    };


    class AreaLight : public Light {
    public:
        AreaLight(const Color& Le, const GeometryPtr& geometry,
            const Transform& toWorld, uint32_t samplesNum);
        ~AreaLight();
        Color sampleL(const Vector3& p, float epsilon, 
            Vector3* wi, Ray* shadowRay) const;
        /*
         * ps: point on area light surface
         * ns: normal on area light surface
         * w: direction L leaving surface
         */
        Color L(const Vector3& ps, const Vector3& ns, const Vector3& w) const;
    private:
        Color mLe;
        Transform mToWorld;
        uint32_t mSamplesNum;
        GeometrySet* mGeometrySet;
    };
}
#endif //GOBLIN_LIGHT_H