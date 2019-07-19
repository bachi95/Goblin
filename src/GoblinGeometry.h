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

        const Vector3& getPosition() const;
        const Vector3& getNormal() const;
        const Vector2& getUV() const;
        const Vector3& getDPDU() const;
        const Vector3& getDPDV() const;
        const Vector3& getDPDX() const;
        const Vector3& getDPDY() const;
        float getDUDX() const;
        float getDVDX() const;
        float getDUDY() const;
        float getDVDY() const;
        Matrix3 getWorldToShade() const;

        void setPosition(const Vector3& position);
        void setNormal(const Vector3& normal);
        void setUV(const Vector2& uv);
        void setDPDU(const Vector3& dpdu);
        void setDPDV(const Vector3& dpdv);
        void setDPDX(const Vector3& dpdx);
        void setDPDY(const Vector3& dpdy);
        void setUVDifferential(float dudx, float dvdx, float dudy, float dvdy);

        void transform(const Transform& t);
    private:
        Vector3 mPosition;
        Vector3 mNormal;
        Vector2 mUV;
        Vector3 mDPDU;
        Vector3 mDPDV;
        Vector3 mDPDX;
        Vector3 mDPDY;
        float mDUDX;
        float mDVDX;
        float mDUDY;
        float mDVDY;
        mutable bool mIsUpdated;
        mutable Matrix3 mWorldToShade;
    };

    inline Fragment::Fragment():
        mDPDX(Vector3::Zero), mDPDY(Vector3::Zero),
        mDUDX(0.0f), mDVDX(0.0f), mDUDY(0.0f), mDVDY(0.0f),
        mIsUpdated(false) {}
    inline const Vector3& Fragment::getPosition() const { return mPosition; }
    inline const Vector3& Fragment::getNormal() const { return mNormal; }
    inline const Vector2& Fragment::getUV() const { return mUV; }
    inline const Vector3& Fragment::getDPDU() const { return mDPDU; }
    inline const Vector3& Fragment::getDPDV() const { return mDPDV; }
    inline const Vector3& Fragment::getDPDX() const { return mDPDX; }
    inline const Vector3& Fragment::getDPDY() const { return mDPDY; }
    inline float Fragment::getDUDX() const { return mDUDX; }
    inline float Fragment::getDVDX() const { return mDVDX; }
    inline float Fragment::getDUDY() const { return mDUDY; }
    inline float Fragment::getDVDY() const { return mDVDY; }

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

    inline void Fragment::setDPDX(const Vector3& dpdx) {
        mIsUpdated = false;
        mDPDX = dpdx;
    }

    inline void Fragment::setDPDY(const Vector3& dpdy) {
        mIsUpdated = false;
        mDPDY = dpdy;
    }

    inline void Fragment::setUVDifferential(float dudx, float dvdx, 
        float dudy, float dvdy) {
        mIsUpdated = false;
        mDUDX = dudx;
        mDVDX = dvdx;
        mDUDY = dudy;
        mDVDY = dvdy;
    }

    struct TriangleIndex {
        unsigned int v[3];
    };

    class Geometry;
    typedef std::shared_ptr<Geometry> GeometryPtr;
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
        for (it = geometryCache.begin(); it != geometryCache.end(); ++it) {
            delete it->second;
        }
        geometryCache.clear();
    }
}

#endif //GOBLIN_GEOMETRY_H
