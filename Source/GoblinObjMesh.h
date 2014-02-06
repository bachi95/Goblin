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
        bool intersect(const Ray& ray) const;
        bool intersect(const Ray& ray, float* epsilon, 
            Fragment* fragment) const;
        float area() const;
        BBox getObjectBound();
        void refine(GeometryList& refinedGeometries);
        bool load();

        bool hasNormal() const;
        bool hasTexCoord() const;
    private:
        void recalculateArea();
    private:
        std::string mFilename;
        BBox mBBox;
        float mArea;
        bool mHasNormal;
        bool mHasTexCoord;
    };

    inline float ObjMesh::area() const { return mArea; }
    inline bool ObjMesh::intersectable() const { return false; }
    inline bool ObjMesh::hasNormal() const { return mHasNormal; }
    inline bool ObjMesh::hasTexCoord() const { return mHasTexCoord; }
}

#endif //GOBLIN_OBJ_MESH
