#include "GoblinDisk.h"
#include "GoblinUtils.h"
#include "GoblinRay.h"
#include "GoblinBBox.h"
#include "GoblinSampler.h"
#include "GoblinScene.h"

namespace Goblin {

Disk::Disk(float r) : mRadius(r) {}

bool Disk::intersect(const Ray& ray, float* epsilon,
	Fragment* fragment) const {
	if (fabs(ray.d.z) < 1e-7f) {
		return false;
	}
	// disk lies on plane z = 0 in local space
	// Oz + Dz * t = 0 -> t = -Oz / Dz;
	float t = -ray.o.z / ray.d.z;
	Vector3 p = ray(t);
	if (t < ray.mint || t > ray.maxt) {
		return false;
	}
	float squareR = p.x * p.x + p.y * p.y;
	if (squareR > mRadius * mRadius) {
		return false;
	}

	ray.maxt = t;
	*epsilon = 1e-3f * t;
	/*
	 * spherical coordinate
	 *
	 * u = Phi / 2PI
	 * v = r / mRadius
	 * x = r * cosPhi = r * cos(2PI * u) = mRadius * v * cosPhi
	 * y = r * sinPhi = r * sin(2PI * u) = mRadius * v * sinPhi
	 * z = 0
	 *
	 * dPx/du = r * -sinPhi * 2PI = -2PI * y
	 * dPy/du = r * cosPhi * 2PI  = 2PI * x
	 * dpdu = [-2PI * y, 2PI * x, 0]
	 *
	 * dPx/dv = mRadius * cosPhi = mRadius * x / r
	 * dPy/dv = mRadius * sinPhi = mRadius * y / r
	 * dpdv = [mRadius * x / r, mRadius * y / r, 0]
	 */
	float r = sqrt(squareR);
	float phi = atan2(p.y, p.x);
	if (phi < 0.0f) {
		phi += TWO_PI;
	}
	float u = phi * INV_TWOPI;
	float v = r / mRadius;
	Vector3 position(p);
	Vector3 normal(Vector3::UnitZ);
	Vector2 uv(u, v);
	Vector3 dpdu(-TWO_PI * p.y, TWO_PI * p.x, 0.0f);
	Vector3 dpdv(mRadius * p.x / r, mRadius * p.y / r, 0.0f);
	*fragment = Fragment(position, normal, uv, dpdu, dpdv);
	return true;
}

bool Disk::occluded(const Ray& ray) const {
	if (fabs(ray.d.z) < 1e-7f) {
		return false;
	}
	// disk lies on plane z = 0 in local space
	// Oz + Dz * t = 0 -> t = -Oz / Dz;
	float t = -ray.o.z / ray.d.z;
	Vector3 p = ray(t);
	if (t < ray.mint || t > ray.maxt) {
		return false;
	}
	return p.x * p.x + p.y * p.y <= mRadius * mRadius;
}

Vector3 Disk::sample(float u1, float u2, Vector3* normal) const {
	*normal = Vector3::UnitZ;
	Vector2 pxy = uniformSampleDisk(u1, u2);
	return Vector3(mRadius * pxy.x, mRadius * pxy.y, 0.0f);
}

BBox Disk::getObjectBound() const {
	return BBox(Vector3(mRadius, mRadius, 0.0f),
		Vector3(-mRadius, -mRadius, 0.0f));
}

Geometry* createDisk(const ParamSet& params, const SceneCache& sceneCache) {
	float radius = params.getFloat("radius", 1.0f);
	return new Disk(radius);
}

} // namespace Goblin