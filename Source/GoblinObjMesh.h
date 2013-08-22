#ifndef GOBLIN_OBJ_MESH
#define GOBLIN_OBJ_MESH

#include "GoblinGeometry.h"
#include "GoblinVertex.h"

#include <string>
#include <vector>

namespace Goblin {
    struct MeshTriangle {
        unsigned int v[3];
    };

    class ObjMesh : public Geometry {
    public:
        ObjMesh(const std::string& filename);
        ~ObjMesh();
        bool init();
        bool load();

        const size_t getVertexNum() const;
        const size_t getFaceNum() const;
        const void* getVertexPtr() const;
        const void* getFacePtr() const;

    public:
        typedef std::vector<Vertex> VertexList;
        typedef std::vector<MeshTriangle> TriangleList;

    private:
        std::string mFilename;
        VertexList mVertices;
        TriangleList mTriangles;
    };

    inline const size_t ObjMesh::getVertexNum() const { 
        return mVertices.size();
    }

    inline const size_t ObjMesh::getFaceNum() const {
        return mTriangles.size();
    }

    inline const void* ObjMesh::getVertexPtr() const {
        return &mVertices[0];
    }

    inline const void* ObjMesh::getFacePtr() const {
        return &mTriangles[0];
    }

}

#endif //GOBLIN_OBJ_MESH_IO