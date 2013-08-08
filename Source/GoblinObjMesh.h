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

    public:
        typedef std::vector<Vertex> VertexList;
        typedef std::vector<MeshTriangle> TriangleList;

        VertexList vertices;
        TriangleList triangles;
    private:
        std::string mFilename;
    };
}

#endif //GOBLIN_OBJ_MESH_IO