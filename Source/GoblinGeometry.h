#ifndef GOBLIN_GEOMETRY_H
#define GOBLIN_GEOMETRY_H
#include "GoblinVertex.h"
#include "GoblinBBox.h"
#include "GoblinUtils.h"

#include <cstdio>
#include <vector>
#include <exception>


namespace Goblin {
    class BBox;
    class Ray;

    //bad.....double reference between Geometry and Primitive
    //TODO change the parameter in Intersect for geometry
    // from Intersection to.... maybe Vertex.....
    class Primitive;
    struct Intersection {
        Intersection(): primitive(NULL) {}

        Vector3 position;
        Vector3 normal;
        const Primitive* primitive;
    };

    struct TriangleIndex {
        unsigned int v[3];
    };

    class Geometry;
    typedef boost::shared_ptr<Geometry> GeometryPtr;
    typedef std::vector<GeometryPtr> GeometryList;

    class Geometry {
    public:
        Geometry();
        virtual ~Geometry() {}
        virtual void init();
        virtual bool intersectable() const;
        virtual bool intersect(const Ray& ray) = 0;
        virtual bool intersect(const Ray& ray, float* epsilon, 
            Intersection* intersection) = 0;
        virtual BBox getObjectBound() = 0;
        virtual void refine(GeometryList& refinedGeometries);
        const size_t getVertexNum() const;
        const size_t getFaceNum() const;
        const void* getVertexPtr(size_t index = 0) const;
        const void* getFacePtr(size_t index = 0) const;
        const size_t getId() const;

    public:
        typedef std::vector<Vertex> VertexList;
        typedef std::vector<TriangleIndex> TriangleList;

    protected:
        static size_t nextGeometryId;
        size_t mGeometryId;

        VertexList mVertices;
        TriangleList mTriangles;
    };

    inline void Geometry::init() {}

    inline bool Geometry::intersectable() const { return true; }

    inline void Geometry::refine(GeometryList& refinedGeometries) {
        std::cerr << "unimplemented Geometry::refine" << std::endl;
        throw std::exception();
    }

    inline const size_t Geometry::getVertexNum() const { 
        return mVertices.size();
    }

    inline const size_t Geometry::getFaceNum() const {
        return mTriangles.size();
    }

    inline const void* Geometry::getVertexPtr(size_t index) const {
        return &mVertices[index];
    }

    inline const void* Geometry::getFacePtr(size_t index) const {
        return &mTriangles[index];
    }
}

#endif //GOBLIN_GEOMETRY_H
