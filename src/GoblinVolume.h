#ifndef GOBLIN_VOLUME_H
#define GOBLIN_VOLUME_H

#include "GoblinBBox.h"
#include "GoblinColor.h"
#include "GoblinFactory.h"
#include "GoblinTransform.h"

namespace Goblin {
class Ray;
class RNG;

class VolumeRegion {
public:
	VolumeRegion(float g, float stepSize, int sampleNum, const BBox& b,
		const Transform& toWorld, bool isHomogeneous) :
		mG(g), mStepSize(stepSize), mSampleNum(sampleNum),
		mLocalRegion(b), mToWorld(toWorld), mIsHomogeneous(isHomogeneous)
	{}

	virtual ~VolumeRegion() {}

	// query volume rendering coefficients at world space position p
	virtual void eval(const Vector3& p, Color& attenuation, Color& scatter,
		Color& emission) const = 0;

	// query extinction coefficient(sigmaT) at world space position p
	// sigmaT = sigmaA + sigmaS
	virtual Color getAttenuation(const Vector3& p) const = 0;

	// calculate transmittance alone input ray
	virtual Color transmittance(const Ray& ray, const RNG& rng) const = 0;

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

	bool isHomogeneous() const { return mIsHomogeneous; }

protected:
	// param for Henyey-Greenstein evaluation
	float mG;
	// sample step size for ray marching calculation
	float mStepSize;
	int mSampleNum;
	BBox mLocalRegion;
	Transform mToWorld;
	bool mIsHomogeneous;
};

// homogeneous volume implementation that absorb/emission/scatter
// coefficients are constant
class HomogeneousVolumeRegion : public VolumeRegion {
public:
	HomogeneousVolumeRegion(const Color& attenuation, const Color& albedo,
		const Color& emission, float g, int sampleNum,
		const BBox& b, const Transform& toWorld):
		VolumeRegion(g, 0.0f, sampleNum, b, toWorld, true),
		mAttenuation(attenuation), mScatter(attenuation * albedo),
		mEmission(emission) {}

	void eval(const Vector3& p, Color& attenuation, Color& scatter,
		Color& emission) const {
		bool inside = mLocalRegion.contain(mToWorld.invertPoint(p));
		if (inside) {
			attenuation = mAttenuation;
			scatter = mScatter;
			emission = mEmission;
		} else {
			attenuation = Color(0.0f);
			scatter = Color(0.0f);
			emission = Color(0.0f);
		}
	}

	Color getAttenuation(const Vector3& p) const {
		bool inside = mLocalRegion.contain(mToWorld.invertPoint(p));
		return inside ? mAttenuation : Color(0.0f);
	}

	Color transmittance(const Ray& ray, const RNG& rng) const;

private:
	// the extinction coefficient (sigmaT)
	Color mAttenuation;
	// the scattering coefficient (sigmaS)
	Color mScatter;
	// the emission coefficient (Lve)
	Color mEmission;
};

class VolumeGrid;

class HeterogeneousVolumeRegion : public VolumeRegion {
public:
	HeterogeneousVolumeRegion(VolumeGrid* density, const Color& albedo, float g,
		float stepSize, int sampleNum, const BBox& b, const Transform& toWorld);

	~HeterogeneousVolumeRegion();

	void eval(const Vector3& p, Color& attenuation, Color& scatter,
		Color& emission) const;

	Color getAttenuation(const Vector3& p) const;

	Color transmittance(const Ray& ray, const RNG& rng) const;

private:
	VolumeGrid* mDensity;
	Color mAlbedo;
};

// Henyey-Greenstein phase function.
inline float phaseHG(const Vector3& wi, const Vector3& wo, float g) {
	if (g < 1e-3) {
		return 0.25f * INV_PI;
	} else {
		float cosTheta = dot(wi, wo);
		return 0.25f * INV_PI * (1.0f - g * g) / 
			powf(1.0f + g * g - 2.0f * g * cosTheta, 1.5f);
	}
}

class ParamSet;
class SceneCache;

class HomogeneousVolumeCreator : public
	Creator<VolumeRegion, const ParamSet&, const SceneCache&> {
public:
	VolumeRegion* create(const ParamSet& params,
		const SceneCache& sceneCache) const;
};

class HeterogeneousVolumeCreator : public
	Creator<VolumeRegion, const ParamSet&, const SceneCache&> {
public:
	VolumeRegion* create(const ParamSet& params,
		const SceneCache& sceneCache) const;
};
}

#endif //GOBLIN_VOLUME_H

