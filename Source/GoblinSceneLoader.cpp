#include "GoblinSceneLoader.h"
#include "GoblinModel.h"
#include "GoblinObjMesh.h"
#include "GoblinSphere.h"
#include "GoblinCamera.h"
#include "GoblinFilm.h"
#include "GoblinUtils.h"

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
    // camera related keywords
    static const char* CAMERA = "camera";
    static const char* FOV = "fov";
    static const char* NEAR_PLANE = "near_plane";
    static const char* FAR_PLANE = "far_plane";
    static const char* FILM = "film";
    static const char* RESOLUTION = "resolution";
    static const char* CROP = "crop";
    // geometry related keywords
    static const char* GEOMETRY = "geometry";
    static const char* NAME = "name";
    static const char* TYPE = "type";
    static const char* MESH = "mesh";
    static const char* SPHERE = "sphere";
    static const char* ORIGIN = "origin";
    static const char* RADIUS = "radius";
    static const char* FILENAME = "file";
    // model related keywords
    static const char* MODEL = "model";
    static const char* POSITION = "position";
    static const char* ORIENTATION = "orientation";
    static const char* SCALE = "scale";
    // material related keywords;
    static const char* MATERIAL = "material";

    static bool getChild(const ptree& pt, const char* key, ptree* child) {
        try {
            *child = pt.get_child(key);
            return true;
        }
        catch (boost::property_tree::ptree_bad_path) {
            return false;
        }
    }

    static float parseFloat(const ptree& pt, const char* key, 
            float fallback = 0.0f) {
        try {
            return pt.get<float>(key);
        }
        catch (boost::property_tree::ptree_bad_path) {
            std::cerr << "value non exist for key " << key << std::endl;
            return fallback;
        }
    }

    static string parseString(const ptree& pt, const char* key,
            const char* fallback = "") {
        try {
            return pt.get<std::string>(key);
        }
        catch (boost::property_tree::ptree_bad_path) {
            std::cerr << "value non exist for key " << key << std::endl;
            return fallback;
        }
    }

    static vector<float> parseFloatArray(const ptree& pt, const char* key) {
        vector<float> rv;
        ptree c;
        if(getChild(pt, key, &c)) {
            for(ptree::const_iterator it = c.begin(); it != c.end(); it++) {
                rv.push_back(lexical_cast<float>(it->second.data()));
            }
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

    static Vector4 parseVector4(const ptree& pt, const char* key) {
        std::vector<float> rv = parseFloatArray(pt, key);
        if(rv.size() != 4) {
            std::cerr << "invalid value for Vector4 " << key << std::endl;
            return Vector4::Zero;
        }
        return Vector4(rv[0], rv[1], rv[2], rv[3]);
    }

    static Quaternion parseQuaternion(const ptree& pt, const char* key) {
        std::vector<float> rv = parseFloatArray(pt, key);
        if(rv.size() != 4) {
            std::cerr << "invalid value for Quaternion " << key << std::endl;
            return Quaternion::Identity;
        }
        return Quaternion(rv[0], rv[1], rv[2], rv[3]);
    }

    Film* parseFilm(const ptree& pt) {
        int xRes = 640;
        int yRes = 480;
        float crop[4] = {0.0f, 1.0f, 0.0f, 1.0f};
        string filename = "ray.png";

        ptree filmTree;
        if(getChild(pt, FILM, &filmTree)) {
            Vector2 res = parseVector2(filmTree, RESOLUTION);
            if(res != Vector2::Zero) {
                xRes = static_cast<int>(res.x);
                yRes = static_cast<int>(res.y);
            }
            Vector4 windowCrop = parseVector4(filmTree, CROP);
            if(windowCrop != Vector4::Zero) {
                for(int i = 0; i < 4; ++i) {
                    crop[i] = windowCrop[i];
                }
            }
            filename = parseString(filmTree, FILENAME, filename.c_str());
        }

        std::cout << "\nfilm" << std::endl;
        std::cout << "-res(" << xRes << ", " << yRes << ")" << std::endl;
        std::cout << "-crop(" << crop[0] << " "<< crop[1] << " "<< 
            crop[2] << " "<< crop[3] << ")" << std::endl;
        std::cout << "-filename: " << filename << std::endl;

        return new Film(xRes, yRes, crop, filename);
    }

    CameraPtr parseCamera(const ptree& pt) {
        Film* film = parseFilm(pt);

        Vector3 position =  Vector3::Zero;
        Quaternion orientation = Quaternion::Identity;
        float zn = 0.0f;
        float zf = 1000.0f;
        float fov = 60.0f;

        ptree cameraTree;
        if(getChild(pt, CAMERA, &cameraTree)) {
            position = parseVector3(cameraTree, POSITION);
            orientation = parseQuaternion(cameraTree, ORIENTATION);
            fov = parseFloat(cameraTree, FOV, fov);
            zn = parseFloat(cameraTree, NEAR_PLANE, zn);
            zf = parseFloat(cameraTree, FAR_PLANE, zf);
        }

        std::cout << "\ncamera" << std::endl;
        std::cout << "-position: " << position << std::endl;
        std::cout << "-orientation: " << orientation << std::endl;
        std::cout << "-fov: " << fov << std::endl;
        std::cout << "-near plane: " << zn << std::endl;
        std::cout << "-far plane: " << zf << std::endl;
        return CameraPtr(new Camera(position, orientation, radians(fov), 
            zn, zf, film));
    }

    void parseGeometry(const ptree& pt, GeometryMap* geometryMap) {
        std::cout <<"\ngeometry" << std::endl;
        string geometryType = parseString(pt, TYPE);

        string name = parseString(pt, NAME);
        GeometryPtr geometry;
        if(geometryType == MESH) {
            string file = parseString(pt, FILENAME);
            geometry = GeometryPtr(new ObjMesh(file));
        } else {
            Vector3 origin = parseVector3(pt, ORIGIN);
            float radius = parseFloat(pt, RADIUS);
            geometry = GeometryPtr(new Sphere(origin, radius));
        }
        geometry->init();
        std::cout << "vertex num: " << geometry->getVertexNum() << std::endl;
        std::cout << "face num: " << geometry->getFaceNum() << std::endl;
        std::pair<string, GeometryPtr> pair(name, geometry);
        geometryMap->insert(pair); 
    }

    ModelPtr parseModel(const ptree& pt, GeometryMap* geometryMap) {
        std::cout << "\nmodel" << std::endl;
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
        try {
            read_json(filename, pt);
        }
        catch(boost::property_tree::json_parser::json_parser_error e) {
            std::cerr <<"error reading scene file " << filename << std::endl;
            std::cerr <<e.what() << std::endl;
            return false;
        }

        GeometryMap geometryMap;

        CameraPtr camera = parseCamera(pt);
        scene->setCamera(camera);

        for(ptree::const_iterator it = pt.begin(); it != pt.end(); it++) {
            std::string key(it->first);
            if(key == MODEL) {
                ModelPtr modelPtr = parseModel(it->second, &geometryMap);
                scene->addModel(modelPtr);
            } else if(key == GEOMETRY) {
                parseGeometry(it->second, &geometryMap);
            }
        }
        
        return true;
    }
}
