#ifndef GOBLIN_OBJ_MESH
#define GOBLIN_OBJ_MESH

#include <string>
#include <vector>
#include "GoblinVertex.h"

namespace Goblin {
    struct MeshTriangle {
        unsigned int v[3];
    };

    class ObjMesh {
    public:
        bool load2(const std::string& filename);
        bool load(const std::string& filename);

    public:
        typedef std::vector<Vertex> VertexList;
        typedef std::vector<MeshTriangle> TriangleList;

        VertexList vertices;
        TriangleList triangles;
    };
}

#endif //GOBLIN_OBJ_MESH_IO