#ifndef GOBLIN_OBJ_MESH
#define GOBLIN_OBJ_MESH

#include "GoblinGeometry.h"
#include "GoblinBBox.h"
#include "GoblinFactory.h"
#include "GoblinTriangle.h"

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
    BBox getObjectBound() const;
    void refine(GeometryList& refinedGeometries) const;
    size_t getVertexNum() const;
    size_t getFaceNum() const;
    const Vertex* getVertexPtr(size_t index) const;
    const TriangleIndex* getFacePtr(size_t index) const;

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
    mutable vector<Triangle> mRefinedMeshes;
    VertexList mVertices;
    TriangleList mTriangles;
};

inline size_t ObjMesh::getVertexNum() const { 
    return mVertices.size();
}

inline size_t ObjMesh::getFaceNum() const {
    return mTriangles.size();
}

inline const Vertex* ObjMesh::getVertexPtr(size_t index) const {
    return &mVertices[index];
}

inline const TriangleIndex* ObjMesh::getFacePtr(size_t index) const {
    return &mTriangles[index];
}

inline float ObjMesh::area() const { return mArea; }

inline bool ObjMesh::intersectable() const { return false; }

inline bool ObjMesh::hasNormal() const { return mHasNormal; }

inline bool ObjMesh::hasTexCoord() const { return mHasTexCoord; }


class ParamSet;
class SceneCache;

class MeshGeometryCreator : 
    public Creator<Geometry, const ParamSet&, const SceneCache&> {
public:
    Geometry* create(const ParamSet& params, 
        const SceneCache& sceneCache) const;
};
}

#endif //GOBLIN_OBJ_MESH
