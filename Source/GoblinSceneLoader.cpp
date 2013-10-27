#include "GoblinSceneLoader.h"
#include "GoblinModel.h"
#include "GoblinObjMesh.h"
#include "GoblinSphere.h"
#include "GoblinCamera.h"
#include "GoblinFilm.h"
#include "GoblinLight.h"
#include "GoblinUtils.h"
#include "GoblinPropertyTree.h"

#include <cassert>
#include <iostream>
#include <string>
#include <map>

namespace Goblin {
    using std::vector;
    using std::string;

    typedef std::map<string, GeometryPtr> GeometryMap;
    typedef std::map<string, PrimitivePtr> PrimitiveMap;
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
    static const char* RADIUS = "radius";
    static const char* FILENAME = "file";
    // model related keywords
    static const char* MODEL = "model";
    // instance related keywords
    static const char* INSTANCE = "instance";
    static const char* POSITION = "position";
    static const char* ORIENTATION = "orientation";
    static const char* SCALE = "scale";
    // light related keywords
    static const char* LIGHT = "light";
    static const char* POINT = "point";
    static const char* DIRECTIONAL = "directional";
    static const char* INTENSITY = "intensity";
    static const char* RADIANCE = "radiance";
    static const char* DIRECTION = "direction";
    // material related keywords;
    static const char* MATERIAL = "material";

    static Vector2 parseVector2(const PropertyTree& pt, const char* key) {
        std::vector<float> rv = pt.parseFloatArray(key);
        if(rv.size() != 2) {
            std::cerr << "invalid value for Vector2 " << key << std::endl;
            return Vector2::Zero;
        }
        return Vector2(rv[0], rv[1]);
    }

    static Vector3 parseVector3(const PropertyTree& pt, const char* key) {
        std::vector<float> rv = pt.parseFloatArray(key);
        if(rv.size() != 3) {
            std::cerr << "invalid value for Vector3 " << key << std::endl;
            return Vector3::Zero;
        }
        return Vector3(rv[0], rv[1], rv[2]);
    }

    static Vector4 parseVector4(const PropertyTree& pt, const char* key) {
        std::vector<float> rv = pt.parseFloatArray(key);
        if(rv.size() != 4) {
            std::cerr << "invalid value for Vector4 " << key << std::endl;
            return Vector4::Zero;
        }
        return Vector4(rv[0], rv[1], rv[2], rv[3]);
    }

    static Color parseColor(const PropertyTree& pt, const char* key) {
        std::vector<float> rv = pt.parseFloatArray(key);
        if(rv.size() != 3) {
            std::cerr << "invalid value for Color " << key << std::endl;
            return Color::White;
        }
        return Color(rv[0], rv[1], rv[2], 1.0f);
    }

    static Quaternion parseQuaternion(const PropertyTree& pt, const char* key) {
        std::vector<float> rv = pt.parseFloatArray(key);
        if(rv.size() != 4) {
            std::cerr << "invalid value for Quaternion " << key << std::endl;
            return Quaternion::Identity;
        }
        return Quaternion(rv[0], rv[1], rv[2], rv[3]);
    }

    Film* parseFilm(const PropertyTree& pt) {
        int xRes = 640;
        int yRes = 480;
        float crop[4] = {0.0f, 1.0f, 0.0f, 1.0f};
        string filename = "ray.png";

        PropertyTree filmTree;
        if(pt.getChild(FILM, &filmTree)) {
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
            filename = filmTree.parseString(FILENAME, filename.c_str());
        }

        std::cout << "\nfilm" << std::endl;
        std::cout << "-res(" << xRes << ", " << yRes << ")" << std::endl;
        std::cout << "-crop(" << crop[0] << " "<< crop[1] << " "<< 
            crop[2] << " "<< crop[3] << ")" << std::endl;
        std::cout << "-filename: " << filename << std::endl;

        return new Film(xRes, yRes, crop, filename);
    }

    CameraPtr parseCamera(const PropertyTree& pt) {
        Film* film = parseFilm(pt);

        Vector3 position =  Vector3::Zero;
        Quaternion orientation = Quaternion::Identity;
        float zn = 0.1f;
        float zf = 1000.0f;
        float fov = 60.0f;

        PropertyTree cameraTree;
        if(pt.getChild(CAMERA, &cameraTree)) {
            position = parseVector3(cameraTree, POSITION);
            orientation = parseQuaternion(cameraTree, ORIENTATION);
            fov = cameraTree.parseFloat(FOV, fov);
            zn = cameraTree.parseFloat(NEAR_PLANE, zn);
            zf = cameraTree.parseFloat(FAR_PLANE, zf);
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

    void parseGeometry(const PropertyTree& pt, GeometryMap* geometryMap) {
        string geometryType = pt.parseString(TYPE);

        string name = pt.parseString(NAME);

        std::cout <<"\ngeometry " << name << std::endl;
        GeometryPtr geometry;
        if(geometryType == MESH) {
            string file = pt.parseString(FILENAME);
            geometry = GeometryPtr(new ObjMesh(file));
        } else {
            float radius = pt.parseFloat(RADIUS);
            geometry = GeometryPtr(new Sphere(radius));
        }
        geometry->init();
        std::cout << "vertex num: " << geometry->getVertexNum() << std::endl;
        std::cout << "face num: " << geometry->getFaceNum() << std::endl;
        std::pair<string, GeometryPtr> pair(name, geometry);
        geometryMap->insert(pair); 

    }

    void parseModel(const PropertyTree& pt, GeometryMap* geometryMap, 
        PrimitiveMap* modelMap) {
        string name = pt.parseString(NAME);

        std::cout << "\nmodel " << name << std::endl;
        string geoName = pt.parseString(GEOMETRY);
        GeometryMap::iterator it = geometryMap->find(geoName);
        if(it == geometryMap->end()) {
            std::cerr << "geometry " << geoName << " not defined!\n";
            return;
        }
        PrimitivePtr model(new Model(it->second));

        if(!model->intersectable()) {
            std::vector<PrimitivePtr> primitives;
            primitives.push_back(model);
            PrimitivePtr aggregate(new Aggregate(primitives));
            std::pair<string, PrimitivePtr> pair(name, aggregate);
            modelMap->insert(pair);
        } else {
            std::pair<string, PrimitivePtr> pair(name, model);
            modelMap->insert(pair);
        }
    }

    void parseInstance(const PropertyTree& pt, PrimitiveMap* modelMap, 
        PrimitiveList* instances) {
        string name = pt.parseString(NAME);
        std::cout << "\ninstance " << name <<std::endl;

        string primitiveName = pt.parseString(MODEL);
        PrimitivePtr primitive;
        PrimitiveMap::iterator it = modelMap->find(primitiveName);
        if(it == modelMap->end()) {
            std::cerr << "model " << primitiveName << " not defined!\n";
            return;
        }

        Vector3 position = parseVector3(pt, POSITION);
        std::cout << "-position: " << position << std::endl;
        Quaternion orientation = parseQuaternion(pt, ORIENTATION);
        std::cout << "-orientation: " << orientation << std::endl;
        Vector3 scale = parseVector3(pt, SCALE);
        std::cout << "-scale: " << scale << std::endl;
        Transform toWorld(position, orientation, scale);

        PrimitivePtr instance(new InstancedPrimitive(toWorld, it->second));
        instances->push_back(instance);
    }

    void parseLight(const PropertyTree& pt, LightList* lights) {
        string lightType = pt.parseString(TYPE);

        string name = pt.parseString(NAME);

        std::cout <<"\nlight " << name << std::endl;
        LightPtr light;
        if(lightType == POINT) {
            Color intensity = parseColor(pt, INTENSITY);
            Vector3 position = parseVector3(pt, POSITION);
            std::cout << "-intensity: " << intensity << std::endl;
            std::cout << "-position: " << position << std::endl;
            light = LightPtr(new PointLight(intensity, position));
        } else if(lightType == DIRECTIONAL) {
            Color radiance = parseColor(pt, RADIANCE);
            Vector3 direction = parseVector3(pt, DIRECTION);
            std::cout << "-radiance: " << radiance << std::endl;
            std::cout << "-direction: " << direction << std::endl;
            light = LightPtr(new DirectionalLight(radiance, direction));
        } else {
            std::cerr << "unrecognized light type " << lightType << std::endl;
        }

        if(light) {
            lights->push_back(light);
        }
    }

    ScenePtr SceneLoader::load(const string& filename) {
        ScenePtr scene;
        PropertyTree pt(filename);

        GeometryMap geometryMap;
        // model name may map to the direct model, or its proxy aggregate
        PrimitiveMap modelMap;
        PrimitiveList instances;
        CameraPtr camera = parseCamera(pt);
        LightList lights;

        const PtreeList& nodes = pt.getChildren(); 
        for(PtreeList::const_iterator it = nodes.begin(); 
            it != nodes.end(); it++) {
            std::string key(it->first);
            if(key == MODEL) {
                parseModel(it->second, &geometryMap, &modelMap);
            } else if(key == GEOMETRY) {
                parseGeometry(it->second, &geometryMap);
            } else if(key == INSTANCE) {
                parseInstance(it->second, &modelMap, &instances);
            } else if(key == LIGHT) {
                parseLight(it->second, &lights);
            }
        }

        PrimitivePtr aggregate(new Aggregate(instances));
        return ScenePtr(new Scene(aggregate, camera, lights));
    }
}
