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
        Renderable(const Matrix4& m, const Geometry* g):
            worldMatrix(m), geometry(g) {}
        Matrix4 worldMatrix;
        const Geometry* geometry;
    };

    typedef std::vector<Renderable> RenderList;

    class Primitive;
    typedef boost::shared_ptr<Primitive> PrimitivePtr;
    typedef vector<const Primitive*> PrimitiveList;

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
        virtual void refine(PrimitiveList& refinedPrimitives) const;
        virtual bool intersect(const Ray& ray) const = 0;
        virtual bool intersect(const Ray& ray, float* epsilon, 
            Intersection* intersection) const = 0;
        virtual BBox getAABB() const = 0;
        virtual const MaterialPtr& getMaterial() const;
        virtual const AreaLight* getAreaLight() const;
        virtual void collectRenderList(RenderList& rList, 
            const Matrix4& m = Matrix4::Identity) const = 0;
        static void clearAllocatedPrimitives();
    protected:
        Primitive() {}
    protected:
        static vector<Primitive*> allocatedPrimitives;
    };


    inline bool Primitive::intersectable() const { return true; }
    inline void Primitive::refine(PrimitiveList& refinedPrimitives) const {
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

    inline void Primitive::clearAllocatedPrimitives() {
        for(size_t i = 0; i < allocatedPrimitives.size(); ++i) {
            delete allocatedPrimitives[i];
            allocatedPrimitives[i] = NULL;
        }
    }

    class InstancedPrimitive : public Primitive{
    public:
        bool intersect(const Ray& ray) const;
        bool intersect(const Ray& ray, float* epsilon, 
            Intersection* intersection) const;
        BBox getAABB() const;
        const Vector3& getPosition() const;
        const Quaternion& getOrientation() const;
        const Vector3& getScale() const;
        const Matrix4& getWorldMatrix();
        void collectRenderList(RenderList& rList, 
            const Matrix4& m = Matrix4::Identity) const;
    private:
        InstancedPrimitive(const Transform& toWorld, 
            const Primitive* primitive);
        InstancedPrimitive(const Vector3& position, 
            const Quaternion& orientation,
            const Vector3& scale, const Primitive* primitive);

    private:
        Transform mToWorld;
        const Primitive* mPrimitive;
    friend class InstancePrimitiveCreator;
    };

    class Aggregate : public Primitive {
    public:
        Aggregate(const PrimitiveList& primitives);
        bool intersect(const Ray& ray) const;
        bool intersect(const Ray& ray, float* epsilon, 
            Intersection* intersection) const;
        void collectRenderList(RenderList& rList, 
            const Matrix4& m = Matrix4::Identity) const;
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
