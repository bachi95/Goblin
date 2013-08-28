#include "GoblinSceneLoader.h"
#include "GoblinModel.h"
#include "GoblinObjMesh.h"
#include "GoblinCamera.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <cassert>
#include <iostream>
#include <string>
#include <map>

namespace Goblin {
    using boost::property_tree::ptree;
    using boost::lexical_cast;
    using std::vector;
    using std::string;

    typedef std::map<string, GeometryPtr> GeometryMap;

    static const char* CAMERA = "camera";
    static const char* RESOLUTION = "resolution";
    static const char* GEOMETRY = "geometry";
    static const char* MODEL = "model";
    static const char* NAME = "name";
    static const char* TYPE = "type";
    static const char* MATERIAL = "material";
    static const char* MESH = "mesh";
    static const char* POSITION = "position";
    static const char* ORIENTATION = "orientation";
    static const char* SCALE = "scale";
    static const char* FILE = "file";

    static string parseString(const ptree& pt, const char* key) {
        try {
            return pt.get<std::string>(key);
        }
        catch (boost::property_tree::ptree_bad_path) {
            std::cerr << "value non exist for key " << key << std::endl;
            return "";
        }
    }

    static vector<float> parseFloatArray(const ptree& pt, const char* key) {
        vector<float> rv;
        const ptree& c = pt.get_child(key);
        for(ptree::const_iterator it = c.begin(); it != c.end(); it++) {
            rv.push_back(lexical_cast<float>(it->second.data()));
        }
        return rv;
    }

    static Vector2 parseVector2(const ptree& pt, const char* key) {
        std::vector<float> rv = parseFloatArray(pt, key);
        if(rv.size() != 2) {
            std::cerr << "invalid value for Vector2 " << key << std::endl;
            return Vector2::Zero;
        }
        return Vector2(rv[0], rv[1]);
    }

    static Vector3 parseVector3(const ptree& pt, const char* key) {
        std::vector<float> rv = parseFloatArray(pt, key);
        if(rv.size() != 3) {
            std::cerr << "invalid value for Vector3 " << key << std::endl;
            return Vector3::Zero;
        }
        return Vector3(rv[0], rv[1], rv[2]);
    }

    static Quaternion parseQuaternion(const ptree& pt, const char* key) {
        std::vector<float> rv = parseFloatArray(pt, key);
        if(rv.size() != 4) {
            std::cerr << "invalid value for Quaternion " << key << std::endl;
            return Quaternion(1, 0, 0, 0);
        }
        return Quaternion(rv[0], rv[1], rv[2], rv[3]);
    }

    CameraPtr parseCamera(const ptree& pt) {
        std::cout << "camera" << std::endl;
        Vector2 filmRes = parseVector2(pt, RESOLUTION);
        std::cout << "-resolution: " << filmRes << std::endl;
        Vector3 position = parseVector3(pt, POSITION);
        std::cout << "-position: " << position << std::endl;
        Quaternion orientation = parseQuaternion(pt, ORIENTATION);
        std::cout << "-orientation: " << orientation << std::endl;
        return CameraPtr(new Camera);
    }

    void parseGeometry(const ptree& pt, GeometryMap* geometryMap) {
        std::cout <<"geometry" << std::endl;
        string geometryType = parseString(pt, TYPE);
        string name = parseString(pt, NAME);
        string file = parseString(pt, FILE);
        GeometryPtr geometry(new ObjMesh(file));
        geometry->init();
        std::cout << "vertex num: " << geometry->getVertexNum() << std::endl;
        std::cout << "face num: " << geometry->getFaceNum() << std::endl;
        std::pair<string, GeometryPtr> pair(name, geometry);
        geometryMap->insert(pair); 
    }

    ModelPtr parseModel(const ptree& pt, GeometryMap* geometryMap) {
        std::cout << "model" << std::endl;
        string geoName = parseString(pt, GEOMETRY);
        GeometryPtr geometry;
        GeometryMap::iterator it = geometryMap->find(geoName);
        if(it == geometryMap->end()) {
            std::cerr << "geometry " << geoName << " not defined!" << std::endl;
        }

        Vector3 position = parseVector3(pt, POSITION);
        std::cout << "-position: " << position << std::endl;
        Quaternion orientation = parseQuaternion(pt, ORIENTATION);
        std::cout << "-orientation: " << orientation << std::endl;
        Vector3 scale = parseVector3(pt, SCALE);
        std::cout << "-scale: " << scale << std::endl;
        Transform toWorld(position, orientation, scale);
        return ModelPtr(new Model(toWorld, it->second));
    }


    bool SceneLoader::load(const string& filename, Scene* scene) {
        ptree pt;
        read_json(filename, pt);
        GeometryMap geometryMap;
        for(ptree::const_iterator it = pt.begin(); it != pt.end(); it++) {
            std::string key(it->first);
            if(key == MODEL) {
                ModelPtr modelPtr = parseModel(it->second, &geometryMap);
                scene->addModel(modelPtr);
            } else if(key == GEOMETRY) {
                parseGeometry(it->second, &geometryMap);
            } else if(key == CAMERA) {
                CameraPtr camera = parseCamera(it->second);
                scene->setCamera(camera);
            } else {
                std::cerr << "unsupport object " << key;
            }
        }
        
        std::cout << "parsing done" << std::endl;
        const ModelList& modelList = scene->getModels();
        for(size_t i = 0; i < modelList.size(); ++i) {
            std::cout << "model " << i << std::endl;
            ModelPtr m =modelList[i];
            std::cout << "position " << m->getPosition() << std::endl;
        }

        return true;
    }
}