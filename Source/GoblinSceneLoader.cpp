#include "GoblinSceneLoader.h"
#include "GoblinRenderer.h"
#include "GoblinModel.h"
#include "GoblinObjMesh.h"
#include "GoblinSphere.h"
#include "GoblinMaterial.h"
#include "GoblinCamera.h"
#include "GoblinFilm.h"
#include "GoblinFilter.h"
#include "GoblinLight.h"
#include "GoblinUtils.h"
#include "GoblinPropertyTree.h"
#include "GoblinBVH.h"
#include "GoblinBBox.h"

#include <cassert>
#include <iostream>
#include <string>
#include <map>

namespace Goblin {
    using std::vector;
    using std::string;

    typedef std::map<string, GeometryPtr> GeometryMap;
    typedef std::map<string, PrimitivePtr> PrimitiveMap;
    typedef std::map<string, MaterialPtr> MaterialMap;
    // camera related keywords
    static const char* CAMERA = "camera";
    static const char* FOV = "fov";
    static const char* NEAR_PLANE = "near_plane";
    static const char* FAR_PLANE = "far_plane";
    static const char* FILM = "film";
    static const char* RESOLUTION = "resolution";
    static const char* CROP = "crop";
    static const char* FILTER = "filter";
    static const char* BOX = "box";
    static const char* TRIANGLE = "triangle";
    static const char* GAUSSIAN = "gaussian";
    static const char* FALLOFF = "falloff";
    static const char* MITCHELL = "mitchell";
    static const char* B = "b";
    static const char* C = "c";
    static const char* FILTER_WIDTH = "width";
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
    static const char* AREA = "area";
    static const char* INTENSITY = "intensity";
    static const char* RADIANCE = "radiance";
    static const char* DIRECTION = "direction";
    // material related keywords;
    static const char* MATERIAL = "material";
    static const char* LAMBERT = "lambert";
    static const char* TRANSPARENT = "transparent";
    static const char* REFLECTION = "Kr";
    static const char* REFRACTION = "Kt";
    static const char* REFRACTION_INDEX = "index";
    static const char* DIFFUSE = "Kd";
    // render setting related keywords;
    static const char* RENDER_SETTING = "render_setting";
    static const char* SAMPLE_PER_PIXEL = "sample_per_pixel";
    static const char* MAX_RAY_DEPTH = "max_ray_depth";
    static const char* SAMPLE_NUM = "sample_num";

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

    static Filter* parseFilter(const PropertyTree& pt) {
        float xWidth = 2;
        float yWidth = 2;
        string filterType = BOX;
        PropertyTree filterTree;
        if(pt.getChild(FILTER, &filterTree)) {
            Vector2 filterWidth = parseVector2(filterTree, FILTER_WIDTH);
            if(filterWidth != Vector2::Zero) {
                xWidth = filterWidth.x;
                yWidth = filterWidth.y;
            }
            filterType = filterTree.parseString(TYPE, BOX);
        }

        std::cout << "\nfilter " << filterType << std::endl;
        std::cout << "-width(" << xWidth << ", " << yWidth << ")" << std::endl;
        Filter* filter;
        if(filterType == BOX) {
            filter = new BoxFilter(xWidth, yWidth);
        } else if(filterType == TRIANGLE) {
            filter = new TriangleFilter(xWidth, yWidth);
        } else if(filterType == GAUSSIAN) {
            float falloff = filterTree.parseFloat(FALLOFF, 2.0f);
            std::cout << "-falloff " << falloff << std::endl;
            filter = new GaussianFilter(xWidth, yWidth, falloff);
        } else if(filterType == MITCHELL) {
            float b = filterTree.parseFloat(B, 1.0f / 3.0f);
            float c = filterTree.parseFloat(C, 1.0f / 3.0f);
            std::cout << "b " << b << " c " << c << std::endl;
            filter = new MitchellFilter(xWidth, yWidth, b, c);
        } else {
            filter = new BoxFilter(xWidth, yWidth);
        }
        return filter;
    }

    static Film* parseFilm(const PropertyTree& pt) {
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
        Filter* filter = parseFilter(pt);

        return new Film(xRes, yRes, crop, filter, filename);
    }

    static CameraPtr parseCamera(const PropertyTree& pt) {
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

    static void parseGeometry(const PropertyTree& pt, GeometryMap* geometryMap) {
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
        BBox bbox = geometry->getObjectBound();
        std::cout << "BBox min: " << bbox.pMin << std::endl;
        std::cout << "BBox max: " << bbox.pMax << std::endl;
        std::pair<string, GeometryPtr> pair(name, geometry);
        geometryMap->insert(pair); 

    }

    static void parseModel(const PropertyTree& pt, GeometryMap* geometryMap, 
        MaterialMap* materialMap, PrimitiveMap* modelMap) {
        string name = pt.parseString(NAME);

        std::cout << "\nmodel " << name << std::endl;
        string geoName = pt.parseString(GEOMETRY);
        GeometryMap::iterator it = geometryMap->find(geoName);
        if(it == geometryMap->end()) {
            std::cerr << "geometry " << geoName << " not defined!\n";
            return;
        }

        string materialName = pt.parseString(MATERIAL);
        MaterialMap::iterator mtlIt = materialMap->find(materialName);
        if(mtlIt == materialMap->end()) {
            std::cerr << "material " << materialName << " not defined!\n";
            return;
        }

        PrimitivePtr model(new Model(it->second, mtlIt->second));

        if(!model->intersectable()) {
            std::vector<PrimitivePtr> primitives;
            primitives.push_back(model);
            PrimitivePtr aggregate(new BVH(primitives, 1, "equal_count"));
            std::pair<string, PrimitivePtr> pair(name, aggregate);
            modelMap->insert(pair);
        } else {
            std::pair<string, PrimitivePtr> pair(name, model);
            modelMap->insert(pair);
        }
    }

    static void parseInstance(const PropertyTree& pt, PrimitiveMap* modelMap, 
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
        BBox bbox = instance->getAABB();
        std::cout << "BBox min: " << bbox.pMin << std::endl;
        std::cout << "BBox max: " << bbox.pMax << std::endl;

        instances->push_back(instance);
    }

    static void parseLight(const PropertyTree& pt, vector<Light*>* lights,
        GeometryMap* geometryMap, PrimitiveList* instances) {
        string lightType = pt.parseString(TYPE);

        string name = pt.parseString(NAME);

        std::cout <<"\nlight " << name << std::endl;
        Light* light;
        if(lightType == POINT) {
            Color intensity = parseColor(pt, INTENSITY);
            Vector3 position = parseVector3(pt, POSITION);
            std::cout << "-intensity: " << intensity << std::endl;
            std::cout << "-position: " << position << std::endl;
            light = new PointLight(intensity, position);
        } else if(lightType == DIRECTIONAL) {
            Color radiance = parseColor(pt, RADIANCE);
            Vector3 direction = parseVector3(pt, DIRECTION);
            std::cout << "-radiance: " << radiance << std::endl;
            std::cout << "-direction: " << direction << std::endl;
            light = new DirectionalLight(radiance, direction);
        } else if(lightType == AREA) {
            Color radiance = parseColor(pt, RADIANCE);
            Vector3 position = parseVector3(pt, POSITION);
            Quaternion orientation = parseQuaternion(pt, ORIENTATION);
            string geoName = pt.parseString(GEOMETRY);
            int sampleNum = pt.parseInt(SAMPLE_NUM);
            std::cout << "-radiance: " << radiance << std::endl;
            std::cout << "-position: " << position << std::endl;
            std::cout << "-orientation: " << orientation << std::endl;
            std::cout << "-geometry: " << geoName << std::endl;
            std::cout << "-samples: " << sampleNum << std::endl;
            GeometryMap::iterator it = geometryMap->find(geoName);
            if(it == geometryMap->end()) {
                std::cerr << "geometry " << geoName << " not defined!\n";
                return;
            }
            // TODO: this cause a problem that we can't run time modify
            // the transform for area light since it's not tied in between
            // instance in scene and the transform in area light itself..
            // need to find a way to improve this part
            Transform toWorld(position, orientation, 
                Vector3(1.0f, 1.0f, 1.0f));

            AreaLight* areaLight = new AreaLight(radiance, it->second, 
                toWorld, sampleNum);
            light = areaLight;
            // and we need to push this geometry into scene so that it can
            // be intersection tested.....this is gonna be messy for the 
            // current awkward parsing mechanics.......
            MaterialPtr mtl(new LambertMaterial(Color::White));
            PrimitivePtr model(new Model(it->second, mtl, areaLight));
            PrimitivePtr instance;
            if(!model->intersectable()) {
                PrimitiveList primitives;
                primitives.push_back(model);
                PrimitivePtr aggregate(new BVH(primitives, 1, "equal_count"));
                instance = 
                    PrimitivePtr(new InstancedPrimitive(toWorld, aggregate));
            } else {
                instance = 
                    PrimitivePtr(new InstancedPrimitive(toWorld, model)); 
            }
            instances->push_back(instance);

        } else {
            std::cerr << "unrecognized light type " << lightType << std::endl;
        }

        if(light) {
            lights->push_back(light);
        }
    }

    static void parseMaterial(const PropertyTree& pt, 
        MaterialMap* materialMap) {
        string materialType = pt.parseString(TYPE);
        string name = pt.parseString(NAME);

        std::cout <<"\n" << materialType <<" material " << name << std::endl;
        MaterialPtr material;
        if(materialType == LAMBERT) {
            Color Kd = parseColor(pt, DIFFUSE);
            std::cout << "-Kd: " << Kd << std::endl;
            material = MaterialPtr(new LambertMaterial(Kd));
        } else if(materialType == TRANSPARENT) {
            Color Kr = parseColor(pt, REFLECTION);
            Color Kt = parseColor(pt, REFRACTION);
            float index = pt.parseFloat(REFRACTION_INDEX, 1.5f);
            std::cout << "-Kr: " << Kr << std::endl;
            std::cout << "-Kt: " << Kt << std::endl;
            std::cout << "refraction index: " << index << std::endl;
            material = MaterialPtr(new TransparentMaterial(Kr, Kt, index));
        } else {
            std::cerr << "undefined material type, use default lambert\n";
            material = MaterialPtr(new LambertMaterial(Color::White));
        }

        std::pair<string, MaterialPtr> pair(name, material);
        materialMap->insert(pair); 
    }

    static void parseRenderSetting(const PropertyTree& pt, 
        RenderSetting* setting) {
        int samplePerPixel = 1;
        int maxRayDepth = 5;

        PropertyTree rt;
        if(pt.getChild(RENDER_SETTING, &rt)) {
            samplePerPixel = rt.parseInt(SAMPLE_PER_PIXEL, samplePerPixel);
            maxRayDepth = rt.parseInt(MAX_RAY_DEPTH, maxRayDepth);
        }
        setting->samplePerPixel = samplePerPixel;
        setting->maxRayDepth = maxRayDepth;
        std::cout << "\nrender setting" << std::endl;
        std::cout << "-sample per pixel " << samplePerPixel << std::endl;
        std::cout << "-max ray depth " << maxRayDepth << std::endl;
    }


    ScenePtr SceneLoader::load(const string& filename, 
        RenderSetting* setting) {
        ScenePtr scene;
        PropertyTree pt;
        if(!pt.read(filename)) {
            return scene;
        }

        GeometryMap geometryMap;
        // model name may map to the direct model, or its proxy aggregate
        PrimitiveMap modelMap;
        MaterialMap materialMap;
        PrimitiveList instances;

        if(setting) {
            parseRenderSetting(pt, setting);
        }
        CameraPtr camera = parseCamera(pt);
        vector<Light*> lights;

        const PtreeList& nodes = pt.getChildren(); 
        for(PtreeList::const_iterator it = nodes.begin(); 
            it != nodes.end(); it++) {
            std::string key(it->first);
            if(key == MODEL) {
                parseModel(it->second, &geometryMap, &materialMap, &modelMap);
            } else if(key == GEOMETRY) {
                parseGeometry(it->second, &geometryMap);
            } else if(key == INSTANCE) {
                parseInstance(it->second, &modelMap, &instances);
            } else if(key == LIGHT) {
                parseLight(it->second, &lights, &geometryMap, &instances);
            } else if(key == MATERIAL) {
                parseMaterial(it->second, &materialMap);
            }
        }

        PrimitivePtr aggregate(new BVH(instances, 1, "equal_count"));
        return ScenePtr(new Scene(aggregate, camera, lights));
    }
}
