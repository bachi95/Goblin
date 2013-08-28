#ifndef GOBLIN_MODEL_H
#define GOBLIN_MODEL_H

#include "GoblinTransform.h"
#include <boost/shared_ptr.hpp>

namespace Goblin {
    class Geometry;
    typedef boost::shared_ptr<Geometry> GeometryPtr;
    class Model {
    public:
        Model(const Transform& toWorld, const GeometryPtr& geometry);
        Model(const Vector3& position, const Quaternion& orientation,
            const Vector3& scale, const GeometryPtr& geometry);
        const Vector3& getPosition() const;
        const Quaternion& getOrientation() const;
        const Vector3& getScale() const;
        const Matrix4& getWorldMatrix();

        const GeometryPtr& getGeometry() const;
    private:
        Transform mToWorld;
        GeometryPtr mGeometry;
        //and Material.....
    };
}

#endif // GOBLIN_MODEL_H