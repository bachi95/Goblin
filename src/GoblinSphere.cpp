#include "GoblinSphere.h"
#include "GoblinUtils.h"
#include "GoblinRay.h"
#include "GoblinBBox.h"
#include "GoblinSampler.h"
#include "GoblinScene.h"

namespace Goblin {

Sphere::Sphere(float r) : mRadius(r) {}

bool Sphere::intersect(const Ray& ray, float* epsilon, 
        Fragment* fragment) const {
    float A = squaredLength(ray.d);
    float B = 2.0f * dot(ray.d, ray.o);
    float C = squaredLength(ray.o) - mRadius * mRadius;
    float tNear, tFar;
    if (!quadratic(A, B, C, &tNear, &tFar)) {
        return false;
    }
    if (tNear > ray.maxt || tFar < ray.mint) {
        return false;
    }
    float tHit = tNear;
    if (tHit < ray.mint) {
        tHit = tFar;
        if (tHit > ray.maxt) {
            return false;
        }
    }
    ray.maxt = tHit;
    *epsilon = 1e-3f * tHit;
    Vector3 pHit = ray(tHit);
    /* 
     * spherical cooridnate
     * x = r * sinTheta * cosPhi
     * y = r * sinTheta * sinPhi
     * z = r * cosTheta
     * 0 < = Phi <= 2PI, 0 <= Theta <=PI
     * u = Phi / 2PI
     * v = Theta / PI
     *
     * dPx/du = d(r * sinTheta * cosPhi)/du = r * sinTheta * d(cos(2PI * u))/du =
     *     r * sinTheta *-sin(2PI * u) * 2PI = -2PI * y
     * dPx/du = d(r * sinTheta * sinPhi)/du = r * sinTheta * d(sin(2PI * u))/du =
     *     r * sinTheta *cos(2PI * u) * 2PI = 2PI * x
     * dPx/dz = d(r * cosTheta)/du = 0
     * dpdu = [-2PI * y, 2PI * x, 0]
     * 
     * dPx/dv = d(r * sinTheta * cosPhi)/dv = r * cosPhi * d(sin(PI * v))/dv =
     *     r * cosPhi * cos(PI * v) * PI = z * PI * cosPhi
     * dPy/dv = d(r * sinTheta * sinPhi)/dv = r * sinPhi * d(sin(PI * v))/dv =
     *    r * sinPhi * cos(PI * v) * PI = z * PI * sinPhi
     * dPz/dv = d(r * cosTheta)/dv = r * d(cos(PI * v))/dv = r * PI * -sinTheta =
     *     -PI * r * sinTheta
     * dpdv = PI * (z * cosPhi, z * sinPhi, -r * sinTheta)
     *
     * cosPhi = x / sqrt(x * x + y * y)
     * sinPhi = y / sqrt(x * x + y * y)
     */
    float phi = atan2(pHit.y, pHit.x); 
    if (phi < 0.0f) {
        phi += TWO_PI;
    }
    float u = phi * INV_TWOPI;
    float theta = acos(pHit.z / mRadius);
    float v = theta * INV_PI;
    float invR = 1.0f / sqrt(pHit.x * pHit.x + pHit.y * pHit.y);
    float cosPhi = pHit.x * invR;
    float sinPhi = pHit.y * invR;

    Vector3 position = pHit;
    Vector3 normal = normalize(pHit);
    Vector2 uv(u, v);
    Vector3 dpdu(-TWO_PI * pHit.y, TWO_PI * pHit.x, 0.0f);
    Vector3 dpdv(PI * Vector3(pHit.z * cosPhi, pHit.z * sinPhi,
        -mRadius * sin(theta)));
    *fragment = Fragment(position, normal, uv, dpdu, dpdv);
    return true;
}

bool Sphere::occluded(const Ray& ray) const {
	float A = squaredLength(ray.d);
	float B = 2.0f * dot(ray.d, ray.o);
	float C = squaredLength(ray.o) - mRadius * mRadius;
	float tNear, tFar;
	if (!quadratic(A, B, C, &tNear, &tFar)) {
		return false;
	}
	if (tNear > ray.maxt || tFar < ray.mint) {
		return false;
	}
	float tHit = tNear;
	if (tHit < ray.mint) {
		tHit = tFar;
		if (tHit > ray.maxt) {
			return false;
		}
	}
	return true;
}

Vector3 Sphere::sample(float u1, float u2, Vector3* normal) const {
    *normal = uniformSampleSphere(u1, u2);
    return mRadius * (*normal);
}

Vector3 Sphere::sample(const Vector3& p, float u1, float u2, 
    Vector3* normal) const {
    float squaredRadius = mRadius * mRadius;
    float squaredDistance = squaredLength(p);
    if (squaredDistance - squaredRadius< 1e-4f) {
        return sample(u1, u2, normal);
    }
    // local space, center is at (0, 0, 0)
    Vector3 zAxis = normalize(-p);
    Vector3 xAxis, yAxis;
    coordinateAxises(zAxis, &xAxis, &yAxis);
    
    float sinThetaMax2 = squaredRadius / squaredDistance;
    float cosThetaMax = sqrt(std::max(0.0f, 1.0f - sinThetaMax2));

    Ray ray(p, 
        uniformSampleCone(u1, u2, cosThetaMax, xAxis, yAxis, zAxis), 1e-3f);
    Fragment fragment; 
    float epsilon;
    Vector3 pHit;
    if (intersect(ray, &epsilon, &fragment)) {
        pHit = fragment.getPosition();
    } else {
        // ray scratches over sphere's surface
        pHit = ray(sqrt(squaredDistance) * cosThetaMax); 
    }
    *normal = normalize(pHit);
    return pHit;
}

float Sphere::pdf(const Vector3& p, const Vector3& wi) const {
    // inside the sphere, return uniform weight
    float squaredDistance = squaredLength(p);
    float squaredRadius = mRadius * mRadius;
    if (squaredDistance - squaredRadius< 1e-4f) {
        return Geometry::pdf(p, wi);
    }
    // outside the sphere, use cone pdf 
    float sinThetaMax2 = squaredRadius / squaredDistance;
    float cosThetaMax = sqrt(std::max(0.0f, 1.0f - sinThetaMax2));
    return uniformConePdf(cosThetaMax);
}

BBox Sphere::getObjectBound() const {
    return BBox(Vector3(mRadius, mRadius, mRadius), 
        Vector3(-mRadius, -mRadius, -mRadius));
}

Geometry* createSphere(const ParamSet& params, const SceneCache& sceneCache) {
    float radius = params.getFloat("radius", 1.0f);
    return new Sphere(radius);
}

}