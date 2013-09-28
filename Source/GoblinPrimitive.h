#ifndef GOBLIN_PRIMITIVE_H
#define GOBLIN_PRIMITIVE_H

#include "GoblinTransform.h"
#include "GoblinGeometry.h"

#include <vector>
#include <boost/shared_ptr.hpp>

namespace Goblin {
    /* temp notes:
    three kinds of primitive: instance, model, aggregate
    model = geometry + material
    instance = transform + any kind of primitive
    aggregate = a collection of primitive

    they all implement:
    intersect
    model::intersect -> return geometry.intersect + material
    instance::intersect -> transform the ray to object space
        then mPrimitive.intersect
    aggregate::intersect -> whatever kind of the space travesel
        ex: kd tree/ BVH or simply naive loop through
    */
    typedef boost::shared_ptr<Geometry> GeometryPtr;

    // temporary workaround for hardware rendering sigh
    // collect a list of renderable from the scene by
    // another temp virtual method collectRenderable
    struct Renderable {
        Renderable(const Matrix4& m, const GeometryPtr& g):
            worldMatrix(m), geometry(g) {}
        Matrix4 worldMatrix;
        GeometryPtr geometry;
    };

    typedef std::vector<Renderable> RenderList;

    class Primitive {
    public:
        bool intersect() { return false; }

        virtual void collectRenderList(RenderList& rList, 
            const Matrix4& m = Matrix4::Identity) = 0;
    };

    typedef boost::shared_ptr<Primitive> PrimitivePtr;
    typedef std::vector<PrimitivePtr> PrimitiveList;

    class InstancedPrimitive : public Primitive{
    public:
        InstancedPrimitive(const Transform& toWorld, 
            const PrimitivePtr& primitive);
        InstancedPrimitive(const Vector3& position, 
            const Quaternion& orientation,
            const Vector3& scale, const PrimitivePtr& primitive);

        const Vector3& getPosition() const;
        const Quaternion& getOrientation() const;
        const Vector3& getScale() const;
        const Matrix4& getWorldMatrix();

        void collectRenderList(RenderList& rList, 
            const Matrix4& m = Matrix4::Identity);
    private:
        Transform mToWorld;
        PrimitivePtr mPrimitive;
    };

    class Aggregate : public Primitive {
    public:
        Aggregate(const PrimitiveList& primitives):
            mPrimitives(primitives) {
        }

        void collectRenderList(RenderList& rList, 
            const Matrix4& m = Matrix4::Identity);
    private:
        PrimitiveList mPrimitives;
    };
}

#endif //GOBLIN_PRIMITIVE_H