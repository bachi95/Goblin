#include "GoblinPropertyTree.h"
#include <boost/property_tree/json_parser.hpp>
#include <boost/lexical_cast.hpp>

namespace Goblin {
    using boost::lexical_cast;

    PropertyTree::PropertyTree(const std::string& filename) {
        try {
            read_json(filename, mPtree);
            for(ptree::const_iterator it = mPtree.begin(); 
                it != mPtree.end(); it++) {
                std::string key(it->first);
                std::pair<std::string, PropertyTree> pair(it->first,
                    PropertyTree(it->second));
                mChildren.push_back(pair);
            }
        }
        catch(boost::property_tree::json_parser::json_parser_error e) {
            std::cerr <<"error reading json file " << filename << std::endl;
            std::cerr <<e.what() << std::endl;
        }
    }

    PropertyTree::PropertyTree(const ptree& pt):
        mPtree(pt) {
        for(ptree::const_iterator it = mPtree.begin(); 
            it != mPtree.end(); it++) {
            std::string key(it->first);
            std::pair<std::string, PropertyTree> pair(it->first,
                PropertyTree(it->second));
            mChildren.push_back(pair);
        }
    }

    const PtreeList& PropertyTree::getChildren() const {
        return mChildren;
    }

    bool PropertyTree::getChild(const char* key, PropertyTree* child) const {
        try {
            ptree c = mPtree.get_child(key);
            *child = PropertyTree(c);
            return true;
        }
        catch (boost::property_tree::ptree_bad_path) {
            return false;
        }
    }

    float PropertyTree::parseFloat(const char* key, float fallback) const {
        try {
            return mPtree.get<float>(key);
        }
        catch (boost::property_tree::ptree_bad_path) {
            std::cerr << "value non exist for key " << key << std::endl;
            return fallback;
        }
    }

    std::string PropertyTree::parseString(const char* key, 
        const char* fallback) const {
        try {
            return mPtree.get<std::string>(key);
        }
        catch (boost::property_tree::ptree_bad_path) {
            std::cerr << "value non exist for key " << key << std::endl;
            return fallback;
        }
    }

    std::vector<float> PropertyTree::parseFloatArray(const char* key) const {
        std::vector<float> rv;
        try {
            ptree c = mPtree.get_child(key);
            for(ptree::const_iterator it = c.begin(); it != c.end(); it++) {
                rv.push_back(lexical_cast<float>(it->second.data()));
            }
        }
        catch (boost::property_tree::ptree_bad_path) {
            std::cerr << "value non exist for key " << key << std::endl;
        }
        return rv;
    }

}
