#include "GoblinLight.h"
#include "GoblinRay.h"

namespace Goblin {
    PointLight::PointLight(const Color& I, const Vector3& P):
    intensity(I), position(P) {
        mParams.setInt("type", Point);
        mParams.setColor("intensity", intensity);
        mParams.setVector3("position", position);
    }

    Color PointLight::Li(const Vector3& p, float epsilon, Vector3* wi,
        Ray* shadowRay) const {
        Vector3 dir = position - p;
        *wi = normalize(dir);
        shadowRay->o = p;
        shadowRay->d = *wi;
        shadowRay->mint = epsilon;
        float squaredDistance = squaredLength(dir);
        shadowRay->maxt = sqrt(squaredDistance);
        return intensity / squaredDistance;
    }

    DirectionalLight::DirectionalLight(const Color& R, const Vector3& D):
    radiance(R), direction(D) {
        mParams.setInt("type", Directional);
        mParams.setColor("radiance", radiance);
        mParams.setVector3("direction", direction);
    }

    Color DirectionalLight::Li(const Vector3& p, float epsilon, 
        Vector3* wi, Ray* shadowRay) const {
        *wi = -direction;
        shadowRay->o = p;
        shadowRay->d = *wi;
        shadowRay->mint = epsilon;
        return radiance;
    }
}