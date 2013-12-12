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
            Fragment* fragment);
        BBox getObjectBound();
        void refine(GeometryList& refinedGeometries);
        bool load();

        bool hasNormal() const;
        bool hasTexCoord() const;
    private:
        std::string mFilename;
        BBox mBBox;
        bool mHasNormal;
        bool mHasTexCoord;
    };

    inline bool ObjMesh::intersectable() const { return false; }
    inline bool ObjMesh::hasNormal() const { return mHasNormal; }
    inline bool ObjMesh::hasTexCoord() const { return mHasTexCoord; }
}

#endif //GOBLIN_OBJ_MESH
