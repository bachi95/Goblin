#include "GoblinPrimitive.h"
#include "GoblinRay.h"

namespace Goblin {

    Color Intersection::Le(const Vector3& outDirection) {
        Vector3 ps = fragment.getPosition();
        Vector3 ns = fragment.getNormal();
        const AreaLight* areaLight = primitive->getAreaLight();
        return areaLight == NULL ? 
            Color::Black : areaLight->L(ps, ns, outDirection);
    }

    const MaterialPtr& Intersection::getMaterial() const {
        return primitive->getMaterial();
    }

    InstancedPrimitive::InstancedPrimitive(const Transform& toWorld, 
        const PrimitivePtr& primitive):
        mToWorld(toWorld), mPrimitive(primitive) {}

    InstancedPrimitive::InstancedPrimitive(const Vector3& position, 
        const Quaternion& orientation,
        const Vector3& scale, const PrimitivePtr& primitive):
        mToWorld(position, orientation, scale), mPrimitive(primitive) {}

    bool InstancedPrimitive::intersect(const Ray& ray) {
        Ray r = mToWorld.invertRay(ray);
        return mPrimitive->intersect(r);
    }

    bool InstancedPrimitive::intersect(const Ray& ray, float* epsilon, 
        Intersection* intersection) {
        Ray r = mToWorld.invertRay(ray);
        bool hit = mPrimitive->intersect(r, epsilon, intersection);
        if(hit) {
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
        const Matrix4& m) {
        mPrimitive->collectRenderList(rList, m * mToWorld.getMatrix());
    }

    Aggregate::Aggregate(const PrimitiveList& primitives):
        mInputPrimitives(primitives) {
        for(size_t i = 0; i < mInputPrimitives.size(); ++i) {
            PrimitivePtr primitive = mInputPrimitives[i];
            if(primitive->intersectable()) {
                mRefinedPrimitives.push_back(mInputPrimitives[i]);
            } else {
                primitive->refine(mRefinedPrimitives);
            }
        }

        for(size_t i = 0; i < mRefinedPrimitives.size(); ++i) {
            mAABB.expand(mRefinedPrimitives[i]->getAABB());
        }
    }
    
    void Aggregate::collectRenderList(RenderList& rList, 
        const Matrix4& m) {
        for(size_t i = 0; i < mInputPrimitives.size(); ++i) {
            mInputPrimitives[i]->collectRenderList(rList, m);
        }
    }

    bool Aggregate::intersect(const Ray& ray) {
        for(size_t i = 0; i < mRefinedPrimitives.size(); ++i) {
            if(mRefinedPrimitives[i]->intersect(ray)) {
                return true;
            }
        }
        return false;
    }

    bool Aggregate::intersect(const Ray& ray, float* epsilon, 
        Intersection* intersection) {
        bool hit = false;
        for(size_t i = 0; i < mRefinedPrimitives.size(); ++i) {
            if(mRefinedPrimitives[i]->intersect(ray, epsilon, 
                intersection)) {
                hit = true;
            }
        }
        return hit;
    }

    BBox Aggregate::getAABB() const {
        return mAABB;
    }
}
