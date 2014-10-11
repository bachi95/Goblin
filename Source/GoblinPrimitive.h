#ifndef GOBLIN_PRIMITIVE_H
#define GOBLIN_PRIMITIVE_H

#include "GoblinTransform.h"
#include "GoblinGeometry.h"
#include "GoblinMaterial.h"
#include "GoblinBBox.h"
#include "GoblinUtils.h"
#include "GoblinLight.h"

#include <vector>
#include <exception>

namespace Goblin {
    class Ray;
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

    class Primitive;
    typedef boost::shared_ptr<Primitive> PrimitivePtr;
    typedef std::vector<PrimitivePtr> PrimitiveList;

    struct Intersection {
        Intersection(): primitive(NULL) {}
        Color Le(const Vector3& outDirection);
        const MaterialPtr& getMaterial() const;
        Fragment fragment;
        const Primitive* primitive;
    };

    class Primitive {
    public:
        virtual ~Primitive() {}
        virtual bool intersectable() const;
        virtual void refine(PrimitiveList& refinedPrimitives);
        virtual bool intersect(const Ray& ray) = 0;
        virtual bool intersect(const Ray& ray, float* epsilon, 
            Intersection* intersection) = 0;
        virtual BBox getAABB() const = 0;
        virtual const MaterialPtr& getMaterial() const;
        virtual const AreaLight* getAreaLight() const;
        virtual void collectRenderList(RenderList& rList, 
            const Matrix4& m = Matrix4::Identity) = 0;
    };


    inline bool Primitive::intersectable() const { return true; }
    inline void Primitive::refine(PrimitiveList& refinedPrimitives) {
        std::cerr << "unimplemented Primitive::refine" << std::endl; 
        throw std::exception();
    }

    // TODO meh......instance should have the right to do 
    //material override, add it in later...
    inline const MaterialPtr& Primitive::getMaterial() const {
        std::cerr << "unimplemented Primitive::getMaterial" << std::endl; 
        throw std::exception();
    }

    inline const AreaLight* Primitive::getAreaLight() const {
        std::cerr << "unimplemented Primitive::getAreaLight" << std::endl; 
        throw std::exception();
    }

    class InstancedPrimitive : public Primitive{
    public:
        InstancedPrimitive(const Transform& toWorld, 
            const PrimitivePtr& primitive);
        InstancedPrimitive(const Vector3& position, 
            const Quaternion& orientation,
            const Vector3& scale, const PrimitivePtr& primitive);
        bool intersect(const Ray& ray);
        bool intersect(const Ray& ray, float* epsilon, 
            Intersection* intersection);
        BBox getAABB() const;
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
        Aggregate(const PrimitiveList& primitives);
        bool intersect(const Ray& ray);
        bool intersect(const Ray& ray, float* epsilon, 
            Intersection* intersection);
        void collectRenderList(RenderList& rList, 
            const Matrix4& m = Matrix4::Identity);
        BBox getAABB() const;
    protected:
        PrimitiveList mInputPrimitives;
        PrimitiveList mRefinedPrimitives;
        BBox mAABB;
    };


    class ParamSet;
    class SceneCache;

    class InstancePrimitiveCreator : public 
        Creator<Primitive , const ParamSet&, const SceneCache&> {
    public:
        Primitive* create(const ParamSet& params, 
            const SceneCache& sceneCache) const;
    };

}

#endif //GOBLIN_PRIMITIVE_H
