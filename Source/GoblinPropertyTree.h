#ifndef GOBLIN_PROPERTY_TREE_H
#define GOBLIN_PROPERTY_TREE_H

#include <string>
#include <vector>
#include <boost/property_tree/ptree.hpp>

namespace Goblin {
    class PropertyTree;
    using boost::property_tree::ptree;
    typedef std::vector<std::pair<std::string, PropertyTree>> PtreeList;

    class PropertyTree {
    public:
        PropertyTree(const std::string& fileName);
        PropertyTree(const ptree& pt);
        PropertyTree() {};
        const PtreeList& getChildren() const;
        bool getChild(const char* key, PropertyTree* child) const;
        float parseFloat(const char* key, float fallback = 0.0f) const;
        std::string parseString(const char* key, const char* fallback = "") const;
        std::vector<float> parseFloatArray(const char* key) const;

    private:
        ptree mPtree;
        PtreeList mChildren;
    };
}

#endif //GOBLIN_PROPERTY_TREE_H