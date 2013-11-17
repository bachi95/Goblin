#include "GoblinBVH.h"
#include "GoblinRay.h"
#include "GoblinUtils.h"
#include <iostream>

namespace Goblin {

    struct BVHPrimitiveInfo {
        BVHPrimitiveInfo(const BBox& b, int i):
            bbox(b), primitiveIndexNum(i), center(0.5f * (b.pMin + b.pMax)) {}
        BBox bbox;
        int primitiveIndexNum;
        Vector3 center;
    };

    struct MidComparator {
        MidComparator(int d, float m): dim(d), midPoint(m) {}
        int dim;
        float midPoint;
        bool operator()(const BVHPrimitiveInfo& b) const {
            return b.center[dim] < midPoint;
        }
    };

    struct PointsComparator {
        PointsComparator(int d): dim(d) {}
        int dim;
        bool operator()(const BVHPrimitiveInfo& a, 
            const BVHPrimitiveInfo& b) const {
            return a.center[dim] < b.center[dim];
        }
    };

    struct BVHTreeNode {
        BVHTreeNode() { children[0] = children[1] = NULL; }
        ~BVHTreeNode() {
            if(children[0] != NULL) {
                delete children[0];
                children[0] = NULL;
            }
            if(children[1] != NULL) {
                delete children[0];
                children[1] = NULL;
            }        
        }

        void initLeaf(int first, int n, const BBox& b) {
            firstPrimIndex = first;
            primitivesNum = n;
            bbox = b;
        }

        void initInterior(int axis, BVHTreeNode* c0, BVHTreeNode* c1) {
            children[0] = c0;
            children[1] = c1;
            bbox = c0->bbox;
            bbox.expand(c1->bbox);
            splitAxis = axis;
            primitivesNum = 0;
        }

        BVHTreeNode* children[2];
        BBox bbox;
        int primitivesNum;
        // for interior node
        int splitAxis;
        // for leaf node
        int firstPrimIndex;
    };

    BVH::BVH(const PrimitiveList& primitives, int maxPrimitivesNum,
        const std::string& splitMethod):
        Aggregate(primitives),
        mMaxPrimitivesNum(maxPrimitivesNum),
        mBVHRoot(NULL) {
        std::cout << "specified split method " << splitMethod << std::endl;
        if(splitMethod == "middle") {
            mSplitMethod = Middle;
        } else if(splitMethod == "equal_count") {
            mSplitMethod = EqualCount;
        } else if(splitMethod == "sah") {
            mSplitMethod = SAH;
            std::cout << "sah split\n";
        } else {
            mSplitMethod = SAH;
        }
        // collect BVHPrimitiveInfo list for the recusive BVH construction
        std::vector<BVHPrimitiveInfo> buildInfoList;
        buildInfoList.reserve(mRefinedPrimitives.size());
        for(size_t i = 0; i < mRefinedPrimitives.size(); ++i) {
            BBox b = mRefinedPrimitives[i]->getAABB();
            buildInfoList.push_back(BVHPrimitiveInfo(b, i));
        } 
        
        //for(size_t i = 0; i < buildInfoList.size(); ++i) {
        //    std::cout << 
        //        buildInfoList[i].primitiveIndexNum <<
        //        " c " <<buildInfoList[i].center <<
        //        " pMin " <<buildInfoList[i].bbox.pMin <<
        //        " pMax " <<buildInfoList[i].bbox.pMax << std::endl;
        //}
        PrimitiveList orderedPrims;
        orderedPrims.reserve(mRefinedPrimitives.size());
        int totalNodes = 0;
        mBVHRoot = buildBVHTree(buildInfoList, 0, buildInfoList.size(),
            &totalNodes, orderedPrims);

        //bool verified = true;
        //for(size_t i = 0; i < buildInfoList.size(); ++i) {
        //    int p = buildInfoList[i].primitiveIndexNum;
        //    BBox o = orderedPrims[i]->getAABB();
        //    BBox r = mRefinedPrimitives[p]->getAABB();
        //    if(o.pMin != r.pMin || o.pMax != r.pMax) {
        //        std::cout << "ERROR!! " << p << std::endl;
        //        verified = false;
        //    }
        //}
        //if(verified) {
        //    std::cout << "ordered prims correctly DFS ordered\n";
        //} else {
        //    std::cout << "ERROR!! ordered prims not correctly ordered\n";
        //}
        mRefinedPrimitives.swap(orderedPrims);
        mBVHNodes.reserve(totalNodes);
        uint32_t offset = 0;
        linearizeBVHTree(mBVHRoot, &offset);
        //compactSummary();
    }

    BVH::~BVH() {
        if(mBVHRoot != NULL) {
            delete mBVHRoot;
            mBVHRoot = NULL;
        }
        std::cout << "BVH deleted\n";
    }

    BVHTreeNode* BVH::buildBVHTree(std::vector<BVHPrimitiveInfo> &buildData,
        int start, int end, int* totalNodes, PrimitiveList& orderedPrims) {
        (*totalNodes)++;
        BVHTreeNode* node = new BVHTreeNode();
        BBox bbox;
        for(int i = start; i < end; ++i) {
            bbox.expand(buildData[i].bbox);
        }
        int primitivesNum = end - start;
        // leaf node case
        if(primitivesNum == 1) {
            int firstPrimIndex = orderedPrims.size();
            //leafSummary(buildData, start, end, firstPrimIndex, primitivesNum);
            for(int i = start; i < end; ++i) {
                int pIndex = buildData[i].primitiveIndexNum;
                orderedPrims.push_back(mRefinedPrimitives[pIndex]);
            }
            node->initLeaf(firstPrimIndex, primitivesNum, bbox);
        } else {
            BBox centersUnion;
            for(int i = start; i < end; ++i) {
                centersUnion.expand(buildData[i].center);
            }
            // pick the axis with largest variant to split
            int dim = centersUnion.longestAxis();
            // all primitives clutter in one point... should be a rare case
            // just make this a leaf node then
            if(centersUnion.pMin[dim] == centersUnion.pMax[dim]) {
                int firstPrimIndex = orderedPrims.size();
                //leafSummary(buildData, start, end, 
                //    firstPrimIndex, primitivesNum);
                for(int i = start; i < end; ++i) {
                    int pIndex = buildData[i].primitiveIndexNum;
                    orderedPrims.push_back(mRefinedPrimitives[pIndex]);
                }
                node->initLeaf(firstPrimIndex, primitivesNum, bbox);

                return node;
            }
            int mid = (start + end) / 2; 
            //split interior node by specified split method
            switch (mSplitMethod) {
            case Middle: {
                float midPoint = 0.5f * (centersUnion.pMin[dim] +
                    centersUnion.pMax[dim]);
                BVHPrimitiveInfo* midPtr = std::partition(&buildData[start], 
                    &buildData[end - 1] + 1, MidComparator(dim, midPoint));
                mid = midPtr - &buildData[0];
                break;
            }
            case EqualCount: {
                mid = (start + end) / 2;
                std::nth_element(&buildData[start], &buildData[mid],
                    &buildData[end - 1] + 1, PointsComparator(dim));
                break;
            }
            case SAH: default:
                //std::cout << "split by SAH\n";
                break;
            }
            //splitSummary(buildData, start, end, mid, dim);
            BVHTreeNode* leftChild = buildBVHTree(buildData, 
                start, mid, totalNodes, orderedPrims);
            BVHTreeNode* rightChild = buildBVHTree(buildData, 
                mid, end, totalNodes, orderedPrims);
            node->initInterior(dim, leftChild, rightChild);
        }
        return node;
    }

    uint32_t BVH::linearizeBVHTree(BVHTreeNode* node, uint32_t* offset) {
        mBVHNodes.push_back(CompactBVHNode());
        uint32_t nodeOffset = (*offset)++;
        CompactBVHNode& compactNode = mBVHNodes[nodeOffset];
        compactNode.bbox = node->bbox;

        if(node->primitivesNum > 0) {
            // leaf node case
            compactNode.firstPrimIndex = node->firstPrimIndex;
            compactNode.primitivesNum = node->primitivesNum;
        } else {
            // interior node case
            compactNode.axis = node->splitAxis;
            compactNode.primitivesNum = 0;
            linearizeBVHTree(node->children[0], offset);
            compactNode.secondChildOffset = 
                linearizeBVHTree(node->children[1], offset);
        }
        return nodeOffset;
    }

    // optimized version bbox/ray intersection test by precomputing
    // invDir and using dirIsNeg indexing to avoid swap tMin/tMax
    // if the ray direction is negative
    static inline bool intersect(const BBox& bbox, const Ray& ray,
        const Vector3& invDir, const uint32_t dirIsNeg[3]) {
        // intersect x tabs
        float tMin = (bbox[dirIsNeg[0]].x - ray.o.x) * invDir.x;
        float tMax = (bbox[1 - dirIsNeg[0]].x - ray.o.x) * invDir.x;
        // intersect y tabs
        float tYMin = (bbox[dirIsNeg[1]].y - ray.o.y) * invDir.y;
        float tYMax = (bbox[1 - dirIsNeg[1]].y - ray.o.y) * invDir.y;
        if(tYMax < tMin || tYMin > tMax) {
            return false;
        }
        if(tYMin > tMin) {
            tMin = tYMin;
        }
        if(tYMax < tMax) {
            tMax = tYMax;
        }
        // intersect ztabs
        float tZMin = (bbox[dirIsNeg[2]].z - ray.o.z) * invDir.z;
        float tZMax = (bbox[1 - dirIsNeg[2]].z - ray.o.z) * invDir.z;
        if(tZMax < tMin || tZMin > tMax) {
            return false;
        }
        if(tZMin > tMin) {
            tMin = tZMin;
        }
        if(tZMax < tMax) {
            tMax = tZMax;
        }

        return (tMin < ray.maxt) && (tMax > ray.mint);
    }

    bool BVH::intersect(const Ray& ray) {
        if(mBVHNodes.size() == 0) {
            return false;
        }
        Vector3 invDir(1.0f / ray.d.x, 1.0f / ray.d.y, 1.0f / ray.d.z);
        uint32_t dirIsNeg[3] = {
            ray.d.x < 0.0f, 
            ray.d.y < 0.0f, 
            ray.d.z < 0.0f};
        uint32_t nodeNum = 0;
        uint32_t todoOffset = 0;
        uint32_t todo[64];
        while(true) {
            CompactBVHNode& node = mBVHNodes[nodeNum];
            if(Goblin::intersect(node.bbox, ray, invDir, dirIsNeg)) {
                if(node.primitivesNum > 0) {
                    for(uint32_t i = 0; i < node.primitivesNum; ++i) {
                        uint32_t index = node.firstPrimIndex + i;
                        if(mRefinedPrimitives[index]->intersect(ray)) {
                            return true;
                        }
                    }
                    if(todoOffset == 0) {
                        break;
                    }
                    nodeNum = todo[--todoOffset];
                } else {
                    if(dirIsNeg[node.axis]) {
                        todo[todoOffset++] = nodeNum + 1;
                        nodeNum = node.secondChildOffset;
                    } else {
                        todo[todoOffset++] = node.secondChildOffset;
                        nodeNum = nodeNum + 1;
                    }
                }
            } else {
                if(todoOffset == 0) {
                    break;
                }
                nodeNum = todo[--todoOffset];
            }
        }
        return false;
    }

    bool BVH::intersect(const Ray& ray, float* epsilon, 
        Intersection* intersection) {
        if(mBVHNodes.size() == 0) {
            return false;
        }
        Vector3 invDir(1.0f / ray.d.x, 1.0f / ray.d.y, 1.0f / ray.d.z);
        uint32_t dirIsNeg[3] = {
            ray.d.x < 0.0f, 
            ray.d.y < 0.0f, 
            ray.d.z < 0.0f};
        uint32_t nodeNum = 0;
        uint32_t todoOffset = 0;
        uint32_t todo[64];
        bool hit = false;
        while(true) {
            CompactBVHNode& node = mBVHNodes[nodeNum];
            if(Goblin::intersect(node.bbox, ray, invDir, dirIsNeg)) {
                if(node.primitivesNum > 0) {
                    for(uint32_t i = 0; i < node.primitivesNum; ++i) {
                        uint32_t index = node.firstPrimIndex + i;
                        if(mRefinedPrimitives[index]->intersect(ray, 
                            epsilon, intersection)) {
                            hit = true;
                        }
                    }
                    if(todoOffset == 0) {
                        break;
                    }
                    nodeNum = todo[--todoOffset];
                } else {
                    if(dirIsNeg[node.axis]) {
                        todo[todoOffset++] = nodeNum + 1;
                        nodeNum = node.secondChildOffset;
                    } else {
                        todo[todoOffset++] = node.secondChildOffset;
                        nodeNum = nodeNum + 1;
                    }
                }
            } else {
                if(todoOffset == 0) {
                    break;
                }
                nodeNum = todo[--todoOffset];
            }
        }
        return hit;
    }

    void BVH::leafSummary(const std::vector<BVHPrimitiveInfo> &buildData,
        int start, int end, int firstPrimIndex, int primitivesNum) {
        std::cout << "--------------------------------\n";
        std::cout << "leaf push in :";
        for(int i = start; i < end; ++i) {
            int pIndex = buildData[i].primitiveIndexNum;
            std::cout << pIndex << " ";
        }
        std::cout << " first " << firstPrimIndex << " num " << 
            primitivesNum << std::endl;
    }

    void BVH::splitSummary(const std::vector<BVHPrimitiveInfo> &buildData, 
        int start, int end, int mid, int dim) {
        std::cout << "--------------------------------\n";
        std::cout << "split summary \n";
        std::cout << "split axes " << (char)('x' + dim) << std::endl;
        std::cout << "start " << start << " end " 
            << end  << " mid " << mid << std::endl; 
        std::cout << "left child: ";
        for(int i = start; i < mid; ++i) {
            std::cout << buildData[i].primitiveIndexNum << " ";
        }
        std::cout << std::endl;

        std::cout << "right child: ";
        for(int i = mid; i < end; ++i) {
            std::cout << buildData[i].primitiveIndexNum << " ";
        }
        std::cout << std::endl;
    }

    void BVH::compactSummary() {
        for(size_t i = 0; i < mRefinedPrimitives.size(); ++i) {
            BBox b = mRefinedPrimitives[i]->getAABB();
            std::cout << i << 
                " c " << 0.5f * (b.pMin + b.pMax) <<
                " pMin " <<b.pMin <<
                " pMax " <<b.pMax << std::endl;
        }

        for(size_t i = 0; i < mBVHNodes.size(); ++ i) {
            CompactBVHNode& node = mBVHNodes[i];
            if(node.primitivesNum == 0) {
                std::cout << i << " axis " << (int)node.axis <<
                    " secondChild " << node.secondChildOffset << std::endl;
            } else {
                std::cout << i << " first " << node.firstPrimIndex <<
                    " primNum " << (int)node.primitivesNum << std::endl;
            }
        }
    }
}
