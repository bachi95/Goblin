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

    BVH::BVH(const PrimitiveList& primitives, int maxPrimitivesNum,
        const std::string& splitMethod):
        Aggregate(primitives),
        mMaxPrimitivesNum(maxPrimitivesNum) {
        if(mRefinedPrimitives.size() == 0) {
            return;
        }
        if(splitMethod == "middle") {
            mSplitMethod = Middle;
        } else if(splitMethod == "equal_count") {
            mSplitMethod = EqualCount;
        } else {
            mSplitMethod = EqualCount;
        }
        // collect BVHPrimitiveInfo list for the recusive BVH construction
        std::vector<BVHPrimitiveInfo> buildInfoList;
        buildInfoList.reserve(mRefinedPrimitives.size());
        for(size_t i = 0; i < mRefinedPrimitives.size(); ++i) {
            BBox b = mRefinedPrimitives[i]->getAABB();
            buildInfoList.push_back(BVHPrimitiveInfo(b, i));
        } 
        //buildDataSummary(buildInfoList);
        PrimitiveList orderedPrims;
        orderedPrims.reserve(mRefinedPrimitives.size());
        mBVHNodes.reserve(2 * mRefinedPrimitives.size() - 1);
        uint32_t offset = 0;
        buildLinearBVH(buildInfoList, 0, buildInfoList.size(),
            &offset, orderedPrims);
        mRefinedPrimitives.swap(orderedPrims);
        //compactSummary();
    }

    BVH::~BVH() {}

    uint32_t BVH::buildLinearBVH(std::vector<BVHPrimitiveInfo> &buildData,
        uint32_t start, uint32_t end, uint32_t* offset, 
        PrimitiveList& orderedPrims) {
        mBVHNodes.push_back(CompactBVHNode());
        uint32_t nodeOffset = (*offset)++;
        CompactBVHNode& node = mBVHNodes[nodeOffset];

        BBox bbox;
        for(uint32_t i = start; i < end; ++i) {
            bbox.expand(buildData[i].bbox);
        }
        uint32_t primitivesNum = end - start;
        // leaf node case
        if(primitivesNum == 1) {
            int firstPrimIndex = orderedPrims.size();
            //leafSummary(buildData, start, end, firstPrimIndex, primitivesNum);
            for(uint32_t i = start; i < end; ++i) {
                uint32_t pIndex = buildData[i].primitiveIndexNum;
                orderedPrims.push_back(mRefinedPrimitives[pIndex]);
            }
            node.initLeaf(bbox, firstPrimIndex, primitivesNum);
        } else {
            BBox centersUnion;
            for(uint32_t i = start; i < end; ++i) {
                centersUnion.expand(buildData[i].center);
            }
            // pick the axis with largest variant to split
            int dim = centersUnion.longestAxis();
            // all primitives clutter in one point... should be a rare case
            // just make this a leaf node then
            if(centersUnion.pMin[dim] == centersUnion.pMax[dim]) {
                uint32_t firstPrimIndex = orderedPrims.size();
                //leafSummary(buildData, start, end, 
                //    firstPrimIndex, primitivesNum);
                for(uint32_t i = start; i < end; ++i) {
                    uint32_t pIndex = buildData[i].primitiveIndexNum;
                    orderedPrims.push_back(mRefinedPrimitives[pIndex]);
                }
                node.initLeaf(bbox, firstPrimIndex, primitivesNum);
                return nodeOffset;
            }
            uint32_t mid = (start + end) / 2; 
            // split interior node by specified split method
            switch (mSplitMethod) {
            case Middle: {
                float midPoint = 0.5f * (centersUnion.pMin[dim] +
                    centersUnion.pMax[dim]);
                BVHPrimitiveInfo* midPtr = std::partition(&buildData[start], 
                    &buildData[end - 1] + 1, MidComparator(dim, midPoint));
                mid = midPtr - &buildData[0];
                // can't split down further with middle method, let the 
                // following split methods handle this case then
                if(start!= mid && end != mid) {
                    break;
                }
            }
            case EqualCount: default: {
                mid = (start + end) / 2;
                std::nth_element(&buildData[start], &buildData[mid],
                    &buildData[end - 1] + 1, PointsComparator(dim));
                break;
            }
            }
            //splitSummary(buildData, start, end, mid, dim);
            buildLinearBVH(buildData, 
                start, mid, offset, orderedPrims);
            uint32_t secondChildOffset = buildLinearBVH(buildData, 
                mid, end, offset, orderedPrims);
            node.initInteror(bbox, secondChildOffset, dim);
        }
        return nodeOffset;
    }

    // optimized version bbox/ray intersection test by precomputing
    // invDir and using dirIsNeg indexing to avoid swap tMin/tMax
    // if the ray direction is negative
    static inline bool intersect(const BBox& bbox, const Ray& ray,
        const Vector3& invDir, uint32_t dirIsNeg[3]) {
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

    bool BVH::intersect(const Ray& ray, IntersectFilter f) const {
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
            const CompactBVHNode& node = mBVHNodes[nodeNum];
            if(Goblin::intersect(node.bbox, ray, invDir, dirIsNeg)) {
                if(node.primitivesNum > 0) {
                    for(uint32_t i = 0; i < node.primitivesNum; ++i) {
                        uint32_t index = node.firstPrimIndex + i;
                        if(mRefinedPrimitives[index]->intersect(ray, f)) {
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
        Intersection* intersection, IntersectFilter f) const {
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
            const CompactBVHNode& node = mBVHNodes[nodeNum];
            if(Goblin::intersect(node.bbox, ray, invDir, dirIsNeg)) {
                if(node.primitivesNum > 0) {
                    for(uint32_t i = 0; i < node.primitivesNum; ++i) {
                        uint32_t index = node.firstPrimIndex + i;
                        if(mRefinedPrimitives[index]->intersect(ray, 
                            epsilon, intersection, f)) {
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


    void BVH::buildDataSummary(
            const std::vector<BVHPrimitiveInfo> &buildData) const {
        std::cout << "--------------------------------\n";
        for(size_t i = 0; i < buildData.size(); ++i) {
            std::cout << 
                buildData[i].primitiveIndexNum <<
                " c " <<buildData[i].center <<
                " pMin " <<buildData[i].bbox.pMin <<
                " pMax " <<buildData[i].bbox.pMax << std::endl;
        }
    }

    void BVH::leafSummary(const std::vector<BVHPrimitiveInfo> &buildData,
        int start, int end, int firstPrimIndex, int primitivesNum) const {
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
        int start, int end, int mid, int dim) const {
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

    void BVH::compactSummary() const {
        std::cout << "--------------------------------\n";
        for(size_t i = 0; i < mRefinedPrimitives.size(); ++i) {
            BBox b = mRefinedPrimitives[i]->getAABB();
            std::cout << i << 
                " c " << 0.5f * (b.pMin + b.pMax) <<
                " pMin " <<b.pMin <<
                " pMax " <<b.pMax << std::endl;
        }

        for(size_t i = 0; i < mBVHNodes.size(); ++ i) {
            const CompactBVHNode& node = mBVHNodes[i];
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
