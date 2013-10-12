#ifndef GOBLIN_OBJ_MESH
#define GOBLIN_OBJ_MESH

#include "GoblinGeometry.h"
#include <string>

namespace Goblin {

    class ObjMesh : public Geometry {
    public:
        ObjMesh(const std::string& filename);
        ~ObjMesh();
        void init();
        bool intersectable() const;
        bool intersect(const Ray& ray);
        void refine(GeometryList& refinedGeometries);
        bool load();

    private:
        std::string mFilename;
    };

    inline bool ObjMesh::intersectable() const { return false; }
}

#endif //GOBLIN_OBJ_MESH
