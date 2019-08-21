#ifndef GOBLIN_POLYGON_MESH
#define GOBLIN_POLYGON_MESH

#include "GoblinGeometry.h"
#include "GoblinBBox.h"
#include "GoblinTriangle.h"

#include <string>

namespace Goblin {

class PolygonMesh : public Geometry {
public:
	PolygonMesh(const std::string& filename);

    ~PolygonMesh() = default;

	bool intersectable() const override {
		return false;
	}

    bool intersect(const Ray& ray, float* epsilon,
		Fragment* fragment) const override;

	bool occluded(const Ray& ray) const override;

	float area() const override {
		return mArea;
	}

    BBox getObjectBound() const override;

    void refine(GeometryList& refinedGeometries) const override;

	const Vertex* getVertexPtr(size_t index) const {
		return &mVertices[index];
	}

	const TriangleIndex* getFacePtr(size_t index) const {
		return &mTriangles[index];
	}

	bool hasNormal() const {
		return mHasNormal;
	}

	bool hasTexCoord() const {
		return mHasTexCoord;
	}

private:
	bool loadObjMesh();

    void recalculateArea();

private:
    std::string mFilename;
    BBox mBBox;
    float mArea;
    bool mHasNormal;
    bool mHasTexCoord;
    std::vector<Triangle> mTriangleGeometries;
    std::vector<Vertex> mVertices;
	std::vector<TriangleIndex> mTriangles;
};

class ParamSet;
class SceneCache;

Geometry* createPolygonMesh(const ParamSet& params,
	const SceneCache& sceneCache);

}

#endif //GOBLIN_POLYGON_MESH
