#include "GoblinParamSet.h"
#include "GoblinPrimitive.h"
#include "GoblinRay.h"
#include "GoblinScene.h"

namespace Goblin {
std::vector<Primitive*> Primitive::allocatedPrimitives;

Color Intersection::Le(const Vector3& outDirection) {
	Vector3 ps = fragment.getPosition();
	Vector3 ns = fragment.getNormal();
	const AreaLight* areaLight = primitive->getAreaLight();
	return areaLight == nullptr ? 
		Color::Black : areaLight->L(ps, ns, outDirection);
}

const Light* Intersection::getLight() const {
	return primitive == nullptr ? nullptr : primitive->getAreaLight();
}

const MaterialPtr& Intersection::getMaterial() const {
	return primitive->getMaterial();
}

bool Intersection::isCameraLens() const {
	return primitive->isCameraLens();
}

bool Intersection::isLight() const {
	return primitive != nullptr && primitive->getAreaLight() != nullptr;
}

void Intersection::computeUVDifferential(const RayDifferential& ray) {
	float dudx, dvdx, dudy, dvdy;
	dudx = dvdx = dudy = dvdy = 0.0f;
	if (ray.hasDifferential) {
		const Vector3& p = fragment.getPosition();
		const Vector3& n = fragment.getNormal();
		/*
			* to get the dx, dy auxiliary ray intersection on plane formed
			* by point p and normal n: ax + by + cz + d = 0, n = (a, b, c)
			* dot(intersection, n) + d = dot(p, n) + d where
			* intersection = Rdx.o + t * Rdx.d =>
			* t = (dot(p, n) - dot(Rdx.o, n)) / dot(Rdx.d, n)
			*/
		float minusD = dot(p, n);
		float tdx = (minusD - dot(ray.dxOrigin, n)) / dot(ray.dxDir, n);
		float tdy = (minusD - dot(ray.dyOrigin, n)) / dot(ray.dyDir, n);
		if (!isNaN(tdx) && !isNaN(tdy)) {
			Vector3 pdx = ray.dxOrigin + tdx * ray.dxDir;
			Vector3 pdy = ray.dyOrigin  +tdy * ray.dyDir;
			Vector3 dpdx = pdx - p;
			Vector3 dpdy = pdy - p;
			fragment.setDPDX(dpdx);
			fragment.setDPDY(dpdy);
			/*
				* solve the equations:
				* dpdx = dpdu * dudx + dpdv * dvdx
				* dpdy = dpdu * dudy + dpdv * dvdy
				* to get dudx, dvdx, dudy, dvdy
				*/
			int axis[2];
			if (fabs(n.x) > fabs(n.y) && fabs(n.x) > fabs(n.z)) {
				// closer to yz plane
				axis[0] = 1;
				axis[1] = 2;
			} else if (fabs(n.y) > fabs(n.z)) {
				// closer to xz plane
				axis[0] = 0;
				axis[1] = 2;
			} else {
				// closer to xy plane
				axis[0] = 0;
				axis[1] = 1;
			}
			const Vector3& dpdu = fragment.getDPDU();
			const Vector3& dpdv = fragment.getDPDV();
			float A[2][2];
			A[0][0] = dpdu[axis[0]];
			A[0][1] = dpdv[axis[0]];
			A[1][0] = dpdu[axis[1]];
			A[1][1] = dpdv[axis[1]];
			float Bx[2];
			Bx[0] = dpdx[axis[0]];
			Bx[1] = dpdx[axis[1]];
			if (!solve2x2LinearSystem(A, Bx, &dudx, &dvdx)) {
				dudx = dvdx = 0.0f;
			}
			float By[2];
			By[0] = dpdy[axis[0]];
			By[1] = dpdy[axis[1]];
			if (!solve2x2LinearSystem(A, By, &dudy, &dvdy)) {
				dudy = dvdy = 0.0f;
			}
		}
	}
	fragment.setUVDifferential(dudx, dvdx, dudy, dvdy);
}

InstancedPrimitive::InstancedPrimitive(const Transform& toWorld, 
	const Primitive* primitive):
	mToWorld(toWorld), mPrimitive(primitive) {}

InstancedPrimitive::InstancedPrimitive(const Vector3& position, 
	const Quaternion& orientation,
	const Vector3& scale, const Primitive* primitive):
	mToWorld(position, orientation, scale), mPrimitive(primitive) {}

bool InstancedPrimitive::intersect(const Ray& ray, 
	IntersectFilter f) const {
	Ray r = mToWorld.invertRay(ray);
	return mPrimitive->intersect(r, f);
}

bool InstancedPrimitive::intersect(const Ray& ray, float* epsilon, 
	Intersection* intersection, IntersectFilter f) const {
	Ray r = mToWorld.invertRay(ray);
	bool hit = mPrimitive->intersect(r, epsilon, intersection, f);
	if (hit) {
		intersection->fragment.transform(mToWorld);
		ray.maxt = r.maxt;
	}
	return hit;
}


BBox InstancedPrimitive::getAABB() const {
	return mToWorld.onBBox(mPrimitive->getAABB());
}

const Vector3& InstancedPrimitive::getPosition() const {
	return mToWorld.getPosition();
}

const Quaternion& InstancedPrimitive::getOrientation() const {
	return mToWorld.getOrientation();
}

const Vector3& InstancedPrimitive::getScale() const {
	return mToWorld.getScale();
}

const Matrix4& InstancedPrimitive::getWorldMatrix() {
	return mToWorld.getMatrix();
}

void InstancedPrimitive::collectRenderList(RenderList& rList, 
	const Matrix4& m) const {
	mPrimitive->collectRenderList(rList, m * mToWorld.getMatrix());
}

Aggregate::Aggregate(const PrimitiveList& primitives):
	mInputPrimitives(primitives) {
	for (size_t i = 0; i < mInputPrimitives.size(); ++i) {
		const Primitive* primitive = mInputPrimitives[i];
		if (primitive->intersectable()) {
			mRefinedPrimitives.push_back(mInputPrimitives[i]);
		} else {
			primitive->refine(mRefinedPrimitives);
		}
	}

	for (size_t i = 0; i < mRefinedPrimitives.size(); ++i) {
		mAABB.expand(mRefinedPrimitives[i]->getAABB());
	}
}
    
void Aggregate::collectRenderList(RenderList& rList, 
	const Matrix4& m) const {
	for (size_t i = 0; i < mInputPrimitives.size(); ++i) {
		mInputPrimitives[i]->collectRenderList(rList, m);
	}
}

bool Aggregate::intersect(const Ray& ray, IntersectFilter f) const {
	for (size_t i = 0; i < mRefinedPrimitives.size(); ++i) {
		if (f != nullptr && !f(mRefinedPrimitives[i], ray)) {
			return false;
		}
		if (mRefinedPrimitives[i]->intersect(ray, f)) {
			return true;
		}
	}
	return false;
}

bool Aggregate::intersect(const Ray& ray, float* epsilon, 
	Intersection* intersection, IntersectFilter f) const {
	bool hit = false;
	for (size_t i = 0; i < mRefinedPrimitives.size(); ++i) {
		if (f != nullptr && !f(mRefinedPrimitives[i], ray)) {
			return false;
		}
		if (mRefinedPrimitives[i]->intersect(ray, epsilon, 
			intersection, f)) {
			hit = true;
		}
	}
	return hit;
}

BBox Aggregate::getAABB() const {
	return mAABB;
}

Primitive* createInstance(const ParamSet& params,
	const SceneCache& sceneCache) {
	std::string primitiveName = params.getString("model");
	const Primitive* primitive = sceneCache.getPrimitive(primitiveName);
	Transform toWorld = getTransform(params);
	Primitive* instance = new InstancedPrimitive(toWorld, primitive);
	Primitive::allocatedPrimitives.push_back(instance);
	return instance;
}

}
