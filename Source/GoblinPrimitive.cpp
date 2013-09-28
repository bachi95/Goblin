#include "GoblinPrimitive.h"

namespace Goblin {
    InstancedPrimitive::InstancedPrimitive(const Transform& toWorld, 
        const PrimitivePtr& primitive):
        mToWorld(toWorld), mPrimitive(primitive) {}

    InstancedPrimitive::InstancedPrimitive(const Vector3& position, 
        const Quaternion& orientation,
        const Vector3& scale, const PrimitivePtr& primitive):
        mToWorld(position, orientation, scale), mPrimitive(primitive) {}

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

    void Aggregate::collectRenderList(RenderList& rList, 
        const Matrix4& m) {
        for(size_t i = 0; i < mPrimitives.size(); ++i) {
            mPrimitives[i]->collectRenderList(rList, m);
        }
    }

}