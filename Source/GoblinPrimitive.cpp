#include "GoblinParamSet.h"
#include "GoblinPrimitive.h"
#include "GoblinRay.h"
#include "GoblinScene.h"

namespace Goblin {
    vector<Primitive*> Primitive::allocatedPrimitives;

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
        const Primitive* primitive):
        mToWorld(toWorld), mPrimitive(primitive) {}

    InstancedPrimitive::InstancedPrimitive(const Vector3& position, 
        const Quaternion& orientation,
        const Vector3& scale, const Primitive* primitive):
        mToWorld(position, orientation, scale), mPrimitive(primitive) {}

    bool InstancedPrimitive::intersect(const Ray& ray) const {
        Ray r = mToWorld.invertRay(ray);
        return mPrimitive->intersect(r);
    }

    bool InstancedPrimitive::intersect(const Ray& ray, float* epsilon, 
        Intersection* intersection) const {
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
        const Matrix4& m) const {
        mPrimitive->collectRenderList(rList, m * mToWorld.getMatrix());
    }

    Aggregate::Aggregate(const PrimitiveList& primitives):
        mInputPrimitives(primitives) {
        for(size_t i = 0; i < mInputPrimitives.size(); ++i) {
            const Primitive* primitive = mInputPrimitives[i];
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
        const Matrix4& m) const {
        for(size_t i = 0; i < mInputPrimitives.size(); ++i) {
            mInputPrimitives[i]->collectRenderList(rList, m);
        }
    }

    bool Aggregate::intersect(const Ray& ray) const {
        for(size_t i = 0; i < mRefinedPrimitives.size(); ++i) {
            if(mRefinedPrimitives[i]->intersect(ray)) {
                return true;
            }
        }
        return false;
    }

    bool Aggregate::intersect(const Ray& ray, float* epsilon, 
        Intersection* intersection) const {
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


    Primitive* InstancePrimitiveCreator::create(const ParamSet& params,
        const SceneCache& sceneCache) const {
        string primitiveName = params.getString("model");
        const Primitive* primitive = sceneCache.getPrimitive(primitiveName);

        Vector3 position = params.getVector3("position");
        Quaternion orientation = getQuaternion(params);
        Vector3 scale = params.getVector3("scale");
        Transform toWorld(position, orientation, scale);
        Primitive* instance = new InstancedPrimitive(toWorld, primitive);
        Primitive::allocatedPrimitives.push_back(instance);
        return instance;
    }

}
