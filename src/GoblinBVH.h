#ifndef GOBLIN_BVH_H
#define GOBLIN_BVH_H
#include "GoblinPrimitive.h"
namespace Goblin {
struct BVHPrimitiveInfo;
struct BVHTreeNode;

struct CompactBVHNode {
    BBox bbox;
    union {
        uint32_t firstPrimIndex; // leaf
        uint32_t secondChildOffset; // interior
    };
    uint8_t primitivesNum;
    uint8_t axis;
    uint8_t pad[2];

    void initLeaf(const BBox& b, uint32_t first, uint8_t n) {
        bbox = b;
        firstPrimIndex = first;
        primitivesNum = n;
    }

    void initInteror(const BBox& b, uint32_t secondChild, uint8_t dim) {
        bbox = b;
        secondChildOffset = secondChild;
        axis = dim;
        primitivesNum = 0;
    }
};

class BVH {
public:
    BVH(const PrimitiveList& primitives, int maxPrimitivesNum,
        const std::string& splitMethod);

	~BVH() = default;

	bool intersect(const Ray& ray, float* epsilon,
		Intersection* intersection, IntersectFilter f) const;

	bool occluded(const Ray& ray, IntersectFilter f) const;

	BBox getAABB() const {
		return mAABB;
	}

private:
    //the BVH we build is a flatten binary tree in DFS order, the node
    //is defined as a compact 32byte class for cache line friendly access
    uint32_t buildLinearBVH(std::vector<BVHPrimitiveInfo> &buildData,
        uint32_t start, uint32_t end, uint32_t* offset,
        PrimitiveList& orderedPrims);

    // these are all just temp debug logging, should find a better verify process
    void buildDataSummary(
        const std::vector<BVHPrimitiveInfo> &buildData) const;

    void leafSummary(const std::vector<BVHPrimitiveInfo> &buildData,
        int start, int end, int firstPrimIndex, int primitivesNum) const ;

    void splitSummary(const std::vector<BVHPrimitiveInfo> &buildData,
        int start, int end, int mid, int dim) const;

    void compactSummary() const;

private:
    enum SplitMethod {
        Middle,
        EqualCount
    };

    int mMaxPrimitivesNum;
    SplitMethod mSplitMethod;
    std::vector<CompactBVHNode> mBVHNodes;
	PrimitiveList mRefinedPrimitives;
	BBox mAABB;
};
}

#endif //GOBLIN_BVH_H