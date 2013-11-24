#ifndef GOBLIN_LIGHT_H
#define GOBLIN_LIGHT_H

#include "GoblinColor.h"
#include "GoblinVector.h"
#include "GoblinParamSet.h"
#include "GoblinUtils.h"

#include <vector>

namespace Goblin {
    class Ray;

    class Light {
    public:
        enum Type {
            Point = 0,
            Directional = 1,
        };
        virtual ~Light() {};
        virtual Color Li(const Vector3& p, float epsilon, 
            Vector3* wi, Ray* shadowRay) const = 0;
        ParamSet& getParams();
    protected:
        ParamSet mParams;
    };

    inline ParamSet& Light::getParams() {
        return mParams;
    }

    typedef boost::shared_ptr<Light> LightPtr;
    typedef std::vector<LightPtr> LightList;

    class PointLight : public Light {
    public:
        PointLight(const Color& intensity, const Vector3& position);
        Color Li(const Vector3& p, float epsilon, Vector3* wi,
            Ray* shadowRay) const;
        Color intensity;
        Vector3 position;
    };

    class DirectionalLight : public Light {
    public:
        DirectionalLight(const Color& radiance, const Vector3& direction);
        Color Li(const Vector3& p, float epsilon, Vector3* wi,
            Ray* shadowRay) const;
        Color radiance;
        Vector3 direction;
    };
}
#endif //GOBLIN_LIGHT_H