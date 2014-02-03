#include "GoblinLight.h"
#include "GoblinRay.h"
#include "GoblinSampler.h"

namespace Goblin {
    PointLight::PointLight(const Color& I, const Vector3& P):
    intensity(I), position(P) {
        mParams.setInt("type", Point);
        mParams.setColor("intensity", intensity);
        mParams.setVector3("position", position);
    }

    Color PointLight::sampleL(const Vector3& p, float epsilon, Vector3* wi,
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

    Color DirectionalLight::sampleL(const Vector3& p, float epsilon, 
        Vector3* wi, Ray* shadowRay) const {
        *wi = -direction;
        shadowRay->o = p;
        shadowRay->d = *wi;
        shadowRay->mint = epsilon;
        return radiance;
    }

    GeometrySet::GeometrySet(const GeometryPtr& geometry):
        mAreaDistribution(NULL) {
        if(geometry->intersectable()) {
            mGeometries.push_back(geometry);
        } else {
            geometry->refine(mGeometries);
        }
        //TODO calculate out the sum area
        //TODO form the area distributioin based on the areas
    }

    GeometrySet::~GeometrySet() {
        if(mAreaDistribution != NULL ) {
            delete mAreaDistribution;
            mAreaDistribution = NULL;
        }
    }

    AreaLight::AreaLight(const Color& Le, const GeometryPtr& geometry,
        const Transform& toWorld, uint32_t samplesNum): mLe(Le), 
        mToWorld(toWorld), mSamplesNum(samplesNum) {
        // TODO geometry init
        mGeometrySet = new GeometrySet(geometry);
    }

    AreaLight::~AreaLight() {
        if(mGeometrySet != NULL) {
            delete mGeometrySet;
            mGeometrySet = NULL;
        }
    }

    Color AreaLight::L(const Vector3& ps, const Vector3& ns, 
        const Vector3& w) const {
        return dot(ns, w) > 0.0f ? mLe : Color::Black;
    }

    Color AreaLight::sampleL(const Vector3& p, float epsilon,
        Vector3* wi, Ray* shadowRay) const {
        Vector3 ns;
        //Vector3 ps = geometries.sample(p, sample, &ns);
        //*wi = normalize(ps - p);
        //*pdf = geometries.pdf(p, *wi);
        //shadowRay->o = p;
        //shadowRay->d = *wi;
        //shadowRay->mint = epsilon;
        //shadowRay->maxt = sqrt(ps - p);
        //return L(ps, ns, -*wi);
        return Color::Black;
    }
}