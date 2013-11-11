#ifndef GOBLIN_OBJ_MESH
#define GOBLIN_OBJ_MESH

#include "GoblinGeometry.h"
#include "GoblinBBox.h"
#include <string>

namespace Goblin {

    class ObjMesh : public Geometry {
    public:
        ObjMesh(const std::string& filename);
        ~ObjMesh();
        void init();
        bool intersectable() const;
        bool intersect(const Ray& ray);
        bool intersect(const Ray& ray, float* epsilon, 
            Intersection* intersection);
        BBox getObjectBound();
        void refine(GeometryList& refinedGeometries);
        bool load();

    private:
        std::string mFilename;
        BBox mBBox;
    };

    inline bool ObjMesh::intersectable() const { return false; }
}

#endif //GOBLIN_OBJ_MESH
