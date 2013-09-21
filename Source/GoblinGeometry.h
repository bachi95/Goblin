#ifndef GOBLIN_GEOMETRY_H
#define GOBLIN_GEOMETRY_H
#include "GoblinVertex.h"
#include <cstdio>
#include <vector>

namespace Goblin {
    struct MeshTriangle {
        unsigned int v[3];
    };

    class Geometry {
    public:
        Geometry();
        virtual ~Geometry() {};
        virtual void init() = 0;
        const size_t getVertexNum() const;
        const size_t getFaceNum() const;
        const void* getVertexPtr() const;
        const void* getFacePtr() const;
        const size_t getId() const;

    public:
        typedef std::vector<Vertex> VertexList;
        typedef std::vector<MeshTriangle> TriangleList;

    protected:
        static size_t nextGeometryId;
        size_t mGeometryId;

        VertexList mVertices;
        TriangleList mTriangles;
    };

    inline const size_t Geometry::getVertexNum() const { 
        return mVertices.size();
    }

    inline const size_t Geometry::getFaceNum() const {
        return mTriangles.size();
    }

    inline const void* Geometry::getVertexPtr() const {
        return &mVertices[0];
    }

    inline const void* Geometry::getFacePtr() const {
        return &mTriangles[0];
    }
}

#endif //GOBLIN_GEOMETRY_H
