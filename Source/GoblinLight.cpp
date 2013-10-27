#include "GoblinLight.h"

namespace Goblin {
    PointLight::PointLight(const Color& I, const Vector3& P):
    intensity(I), position(P) {
        mParams.setInt("type", Point);
        mParams.setColor("intensity", intensity);
        mParams.setVector3("position", position);
    }

    Color PointLight::Li(const Vector3& p, float* epsilon, Vector3* wi) const {
        Vector3 dir = position - p;
        *wi = normalize(dir);
        return intensity / squaredLength(dir);
    }

    DirectionalLight::DirectionalLight(const Color& R, const Vector3& D):
    radiance(R), direction(D) {
        mParams.setInt("type", Directional);
        mParams.setColor("radiance", radiance);
        mParams.setVector3("direction", direction);
    }

    Color DirectionalLight::Li(const Vector3& p, float* epsilon, 
        Vector3* wi) const {
        *wi = -direction;
        return radiance;
    }
}