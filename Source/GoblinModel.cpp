#include "GoblinModel.h"

namespace Goblin {
    Model::Model(const Transform& toWorld, const GeometryPtr& geometry):
        mToWorld(toWorld), mGeometry(geometry) {}

    Model::Model(const Vector3& position, const Quaternion& orientation,
        const Vector3& scale, const GeometryPtr& geometry):
        mToWorld(position, orientation, scale), mGeometry(geometry) {}

    const Vector3& Model::getPosition() const {
        return mToWorld.getPosition();
    }

    const Quaternion& Model::getOrientation() const {
        return mToWorld.getOrientation();
    }

    const Vector3& Model::getScale() const {
        return mToWorld.getScale();
    }

    const Matrix4& Model::getWorldMatrix() {
        return mToWorld.getMatrix();
    }

    const GeometryPtr& Model::getGeometry() const {
        return mGeometry;
    }
}