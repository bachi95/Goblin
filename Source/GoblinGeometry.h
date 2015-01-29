#ifndef GOBLIN_GEOMETRY_H
#define GOBLIN_GEOMETRY_H
#include "GoblinVertex.h"
#include "GoblinBBox.h"
#include "GoblinUtils.h"
#include "GoblinMatrix.h"

#include <cstdio>
#include <vector>
#include <exception>


namespace Goblin {
    class BBox;
    class Ray;
    class Transform;

    class Fragment {
    public:
        Fragment();
        Fragment(const Vector3& p, const Vector3& n, const Vector2& uv,
            const Vector3& dpdu, const Vector3& dpdv);

        Vector3 getPosition() const;
        Vector3 getNormal() const;
        Vector3 getDPDU() const;
        Vector3 getDPDV() const;
        Vector2 getUV() const;
        Matrix3 getWorldToShade() const;

        void setPosition(const Vector3& position);
        void setNormal(const Vector3& normal);
        void setDPDU(const Vector3& dpdu);
        void setDPDV(const Vector3& dpdv);
        void setUV(const Vector2& uv);

        void transform(const Transform& t);
    private:
        Vector3 mPosition;
        Vector3 mNormal;
        Vector2 mUV;
        Vector3 mDPDU;
        Vector3 mDPDV;
        mutable bool mIsUpdated;
        mutable Matrix3 mWorldToShade;
    };

    inline Fragment::Fragment():mIsUpdated(false) {}
    inline Vector3 Fragment::getPosition() const { return mPosition; }
    inline Vector3 Fragment::getNormal() const { return mNormal; }
    inline Vector2 Fragment::getUV() const { return mUV; }
    inline Vector3 Fragment::getDPDU() const { return mDPDU; }
    inline Vector3 Fragment::getDPDV() const { return mDPDV; }

    inline void Fragment::setPosition(const Vector3& position) { 
        mIsUpdated = false;
        mPosition = position;
    }

    inline void Fragment::setNormal(const Vector3& normal) {
        mIsUpdated = false;
        mNormal = normal;
    }

    inline void Fragment::setUV(const Vector2& uv) {
        mIsUpdated = false;
        mUV = uv;
    }

    inline void Fragment::setDPDU(const Vector3& dpdu) {
        mIsUpdated = false;
        mDPDU = dpdu;
    }

    inline void Fragment::setDPDV(const Vector3& dpdv) {
        mIsUpdated = false;
        mDPDV = dpdv;
    }


    struct TriangleIndex {
        unsigned int v[3];
    };

    class Geometry;
    typedef boost::shared_ptr<Geometry> GeometryPtr;
    typedef std::vector<const Geometry*> GeometryList;

    class Geometry {
    public:
        Geometry();
        virtual ~Geometry() {}
        virtual void init();
        virtual bool intersectable() const;
        virtual bool intersect(const Ray& ray) const = 0;
        virtual bool intersect(const Ray& ray, float* epsilon, 
            Fragment* fragment) const = 0;
        virtual Vector3 sample(float u1, float u2, Vector3* normal) const;
        virtual Vector3 sample(const Vector3& p, float u1, float u2, 
            Vector3* normal) const;
        // pdf with respect to solid angle from p and extanded by this geometry
        virtual float pdf(const Vector3& p, const Vector3& wi) const; 
        virtual float area() const = 0;
        virtual BBox getObjectBound() const = 0;
        virtual void refine(GeometryList& refinedGeometries) const;

        virtual size_t getVertexNum() const = 0;
        virtual size_t getFaceNum() const = 0;
        virtual const Vertex* getVertexPtr(size_t index = 0) const = 0;
        virtual const TriangleIndex* getFacePtr(size_t index = 0) const = 0;
        const size_t getId() const;

        static void clearGeometryCache();

    public:
        typedef std::vector<Vertex> VertexList;
        typedef std::vector<TriangleIndex> TriangleList;

    protected:
        static size_t nextGeometryId;
        static std::map<size_t, Geometry*> geometryCache;
        size_t mGeometryId;
    };

    inline void Geometry::init() {}

    inline bool Geometry::intersectable() const { return true; }

    inline Vector3 Geometry::sample(float u1, float u2, 
        Vector3* normal) const {
        std::cerr << "unimplemented Geometry::sample(float, float, Vector3*)"
            << std::endl;
        throw std::exception();
    }

    inline void Geometry::refine(GeometryList& refinedGeometries) const {
        std::cerr << "unimplemented Geometry::refine" << std::endl;
        throw std::exception();
    }

    inline void Geometry::clearGeometryCache() {
        std::map<size_t, Geometry*>::iterator it;
        for(it = geometryCache.begin(); it != geometryCache.end(); ++it) {
            delete it->second;
        }
        geometryCache.clear();
    }
}

#endif //GOBLIN_GEOMETRY_H
