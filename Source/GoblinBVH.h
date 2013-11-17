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
    };

    class BVH : public Aggregate {
    public: 
        BVH(const PrimitiveList& primitives, int maxPrimitivesNum = 1,
            const std::string& splitMethod = "middle");
        ~BVH();
        bool intersect(const Ray& ray);
        bool intersect(const Ray& ray, float* epsilon, 
            Intersection* intersection);
    private:
        BVHTreeNode* buildBVHTree(std::vector<BVHPrimitiveInfo> &buildData,
            int start, int end, int* totalNodes, PrimitiveList& orderedPrims);
        uint32_t linearizeBVHTree(BVHTreeNode* node, uint32_t* offset);

        void leafSummary(const std::vector<BVHPrimitiveInfo> &buildData,
            int start, int end, int firstPrimIndex, int primitivesNum);
        void splitSummary(const std::vector<BVHPrimitiveInfo> &buildData, 
            int start, int end, int mid, int dim);
        void compactSummary();

    private:
        enum SplitMethod {
            Middle, 
            EqualCount,
            SAH
        };
        BVHTreeNode* mBVHRoot;
        int mMaxPrimitivesNum;
        SplitMethod mSplitMethod;

        std::vector<CompactBVHNode> mBVHNodes;
    };
}

#endif //GOBLIN_BVH_H
