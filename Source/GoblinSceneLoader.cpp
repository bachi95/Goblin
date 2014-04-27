#include "GoblinSceneLoader.h"
#include "GoblinBBox.h"
#include "GoblinBVH.h"
#include "GoblinCamera.h"
#include "GoblinFilm.h"
#include "GoblinFilter.h"
#include "GoblinLight.h"
#include "GoblinMaterial.h"
#include "GoblinModel.h"
#include "GoblinObjMesh.h"
#include "GoblinPropertyTree.h"
#include "GoblinRenderer.h"
#include "GoblinSphere.h"
#include "GoblinTexture.h"
#include "GoblinUtils.h"

#include <cassert>
#include <iostream>
#include <map>
#include <boost/filesystem.hpp>

namespace Goblin {
    using boost::filesystem::path;
    typedef std::map<string, GeometryPtr> GeometryMap;
    typedef std::map<string, PrimitivePtr> PrimitiveMap;
    typedef std::map<string, MaterialPtr> MaterialMap;
    typedef std::map<string, TexturePtr> TextureMap;

    // camera related keywords
    static const char* CAMERA = "camera";
    static const char* FOV = "fov";
    static const char* NEAR_PLANE = "near_plane";
    static const char* FAR_PLANE = "far_plane";
    static const char* LENS_RADIUS = "lens_radius";
    static const char* FOCAL_DISTANCE = "focal_distance";
    static const char* FILM = "film";
    static const char* RESOLUTION = "resolution";
    static const char* CROP = "crop";
    // filter related keywords
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
    // texture related keywords:
    static const char* TEXTURE = "texture";
    static const char* CONSTANT = "constant";
    static const char* COLOR = "color";
    static const char* IMAGE = "image";
    static const char* MAPPING = "mapping";
    static const char* UV = "uv";
    static const char* OFFSET = "offset";
    static const char* GAMMA = "gamma";
    static const char* ADDRESS = "address";
    static const char* REPEAT = "repeat";
    static const char* CLAMP = "clamp";
    static const char* BORDER = "border";
    // render setting related keywords;
    static const char* RENDER_SETTING = "render_setting";
    static const char* SAMPLE_PER_PIXEL = "sample_per_pixel";
    static const char* MAX_RAY_DEPTH = "max_ray_depth";
    static const char* SAMPLE_NUM = "sample_num";


    class SceneCache {
    public:
        SceneCache(const path& sceneRoot);
        void addGeometry(const string& name, GeometryPtr g);
        void addPrimitive(const string& name, PrimitivePtr p);
        void addMaterial(const string& name, MaterialPtr m);
        void addTexture(const string& name, TexturePtr t);
        void addInstance(const PrimitivePtr& i);
        void addLight(Light* l);
        const GeometryPtr& getGeometry(const string& name) const;
        const PrimitivePtr& getPrimitive(const string& name) const;
        const MaterialPtr& getMaterial(const string& name) const;
        const TexturePtr& getTexture(const string& name) const;
        const PrimitiveList& getInstances() const;
        const vector<Light*>& getLights() const;
        string resolvePath(const string& filename) const;

    private:
        void initDefault();
        GeometryMap mGeometryMap;
        PrimitiveMap mPrimitiveMap;
        MaterialMap mMaterialMap;
        TextureMap mTextureMap;
        PrimitiveList mInstances;
        vector<Light*> mLights;
        path mSceneRoot;
        string mErrorCode;
    };

    SceneCache::SceneCache(const path& sceneRoot): 
        mSceneRoot(sceneRoot),
        mErrorCode("error") {
        initDefault();
    }

    void SceneCache::initDefault() {
        Color errorColor = Color::Magenta;
        TexturePtr errorTexture(new ConstantTexture(errorColor));
        addTexture(mErrorCode, errorTexture);
        MaterialPtr errorMaterial(new LambertMaterial(errorTexture));
        addMaterial(mErrorCode, errorMaterial);
        GeometryPtr errorGeometry(new Sphere(1.0f));
        addGeometry(mErrorCode, errorGeometry);
        PrimitivePtr errorPrimitive(new Model(errorGeometry, errorMaterial));
        addPrimitive(mErrorCode, errorPrimitive);
    }

    void SceneCache::addGeometry(const string& name, GeometryPtr g) {
        std::pair<string, GeometryPtr> pair(name, g);
        mGeometryMap.insert(pair); 
    }

    void SceneCache::addPrimitive(const string& name, PrimitivePtr p) {
        std::pair<string, PrimitivePtr> pair(name, p);
        mPrimitiveMap.insert(pair); 
    }

    void SceneCache::addMaterial(const string& name, MaterialPtr m) {
        std::pair<string, MaterialPtr> pair(name, m);
        mMaterialMap.insert(pair); 
    }

    void SceneCache::addTexture(const string& name, TexturePtr t) {
        std::pair<string, TexturePtr> pair(name, t);
        mTextureMap.insert(pair); 
    }

    void SceneCache::addInstance(const PrimitivePtr& i) {
        mInstances.push_back(i);
    }

    void SceneCache::addLight(Light* l) {
        mLights.push_back(l);
    }

    const GeometryPtr& SceneCache::getGeometry(const string& name) const {
        GeometryMap::const_iterator it = mGeometryMap.find(name);
        if(it == mGeometryMap.end()) {
            std::cerr << "Geometry " << name << " not defined!\n";
            return mGeometryMap.find(mErrorCode)->second;
        }
        return it->second;
    }

    const PrimitivePtr& SceneCache::getPrimitive(const string& name) const {
        PrimitiveMap::const_iterator it = mPrimitiveMap.find(name);
        if(it == mPrimitiveMap.end()) {
            std::cerr << "Primitive " << name << " not defined!\n";
            return mPrimitiveMap.find(mErrorCode)->second;
        }
        return it->second;
    }

    const MaterialPtr& SceneCache::getMaterial(const string& name) const {
        MaterialMap::const_iterator it = mMaterialMap.find(name);
        if(it == mMaterialMap.end()) {
            std::cerr << "Material " << name << " not defined!\n";
            return mMaterialMap.find(mErrorCode)->second;
        }
        return it->second;
    }

    const TexturePtr& SceneCache::getTexture(const string& name) const {
        TextureMap::const_iterator it = mTextureMap.find(name);
        if(it == mTextureMap.end()) {
            std::cerr << "Texture " << name << " not defined!\n";
            return mTextureMap.find(mErrorCode)->second;
        }
        return it->second;
    }

    const PrimitiveList& SceneCache::getInstances() const {
        return mInstances;
    }

    const vector<Light*>& SceneCache::getLights() const {
        return mLights;
    }

    string SceneCache::resolvePath(const string& filename) const {
        path filePath(filename);
        if(filePath.is_absolute()) {
            return filePath.generic_string();
        } else {
            return (mSceneRoot / filename).generic_string();
        }
    }

    static Vector2 parseVector2(const PropertyTree& pt, const char* key,
        const Vector2& fallback = Vector2::Zero) {
        std::vector<float> rv = pt.parseFloatArray(key);
        if(rv.size() != 2) {
            std::cerr << "invalid value for Vector2 " << key << std::endl;
            return fallback;
        }
        return Vector2(rv[0], rv[1]);
    }

    static Vector3 parseVector3(const PropertyTree& pt, const char* key,
        const Vector3& fallback = Vector3::Zero) {
        std::vector<float> rv = pt.parseFloatArray(key);
        if(rv.size() != 3) {
            std::cerr << "invalid value for Vector3 " << key << std::endl;
            return fallback;
        }
        return Vector3(rv[0], rv[1], rv[2]);
    }

    static Vector4 parseVector4(const PropertyTree& pt, const char* key,
        const Vector4& fallback = Vector4::Zero) {
        std::vector<float> rv = pt.parseFloatArray(key);
        if(rv.size() != 4) {
            std::cerr << "invalid value for Vector4 " << key << std::endl;
            return fallback;
        }
        return Vector4(rv[0], rv[1], rv[2], rv[3]);
    }

    static Color parseColor(const PropertyTree& pt, const char* key,
        const Color& fallback = Color::White) {
        std::vector<float> rv = pt.parseFloatArray(key);
        if(rv.size() != 3) {
            std::cerr << "invalid value for Color " << key << std::endl;
            return fallback;
        }
        return Color(rv[0], rv[1], rv[2], 1.0f);
    }

    static Quaternion parseQuaternion(const PropertyTree& pt, const char* key,
        const Quaternion& fallback = Quaternion::Identity) {
        std::vector<float> rv = pt.parseFloatArray(key);
        if(rv.size() != 4) {
            std::cerr << "invalid value for Quaternion " << key << std::endl;
            return fallback;
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

    static Film* parseFilm(const PropertyTree& pt, SceneCache* sceneCache) {
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
        string filePath = sceneCache->resolvePath(filename);

        std::cout << "\nfilm" << std::endl;
        std::cout << "-res(" << xRes << ", " << yRes << ")" << std::endl;
        std::cout << "-crop(" << crop[0] << " "<< crop[1] << " "<< 
            crop[2] << " "<< crop[3] << ")" << std::endl;
        std::cout << "-filepath: " << filePath << std::endl;
        Filter* filter = parseFilter(pt);

        return new Film(xRes, yRes, crop, filter, filePath);
    }

    static CameraPtr parseCamera(const PropertyTree& pt, Film* film) {
        Vector3 position =  Vector3::Zero;
        Quaternion orientation = Quaternion::Identity;
        float zn = 0.1f;
        float zf = 1000.0f;
        float fov = 60.0f;
        float lensRadius = 0.0f;
        float focalDistance = 1e+10f;
        PropertyTree cameraTree;
        if(pt.getChild(CAMERA, &cameraTree)) {
            position = parseVector3(cameraTree, POSITION);
            orientation = parseQuaternion(cameraTree, ORIENTATION);
            fov = cameraTree.parseFloat(FOV, fov);
            zn = cameraTree.parseFloat(NEAR_PLANE, zn);
            zf = cameraTree.parseFloat(FAR_PLANE, zf);
            lensRadius = cameraTree.parseFloat(LENS_RADIUS, lensRadius);
            focalDistance = cameraTree.parseFloat(FOCAL_DISTANCE, 
                focalDistance);
        }

        std::cout << "\ncamera" << std::endl;
        std::cout << "-position: " << position << std::endl;
        std::cout << "-orientation: " << orientation << std::endl;
        std::cout << "-fov: " << fov << std::endl;
        std::cout << "-near plane: " << zn << std::endl;
        std::cout << "-far plane: " << zf << std::endl;
        std::cout << "-lens radius: " << lensRadius << std::endl;
        std::cout << "-focal distance: " << focalDistance << std::endl;
        return CameraPtr(new Camera(position, orientation, radians(fov), 
            zn, zf, lensRadius, focalDistance, film));
    }

    static void parseGeometry(const PropertyTree& pt, SceneCache* sceneCache) {
        string geometryType = pt.parseString(TYPE);
        string name = pt.parseString(NAME);
        std::cout <<"\ngeometry " << name << std::endl;
        GeometryPtr geometry;
        if(geometryType == MESH) {
            string filename = pt.parseString(FILENAME);
            string filePath = sceneCache->resolvePath(filename);
            geometry = GeometryPtr(new ObjMesh(filePath));
        } else if(geometryType == SPHERE){
            float radius = pt.parseFloat(RADIUS);
            geometry = GeometryPtr(new Sphere(radius));
        } else {
            std::cerr << "unknowened geometry type " << geometryType<< "\n";
            return;
        }
        geometry->init();
        std::cout << "vertex num: " << geometry->getVertexNum() << std::endl;
        std::cout << "face num: " << geometry->getFaceNum() << std::endl;
        BBox bbox = geometry->getObjectBound();
        std::cout << "BBox min: " << bbox.pMin << std::endl;
        std::cout << "BBox max: " << bbox.pMax << std::endl;
        sceneCache->addGeometry(name, geometry);
    }

    static void parseModel(const PropertyTree& pt, SceneCache* sceneCache) {
        string name = pt.parseString(NAME);

        std::cout << "\nmodel " << name << std::endl;
        string geoName = pt.parseString(GEOMETRY);
        GeometryPtr geometry = sceneCache->getGeometry(geoName);

        string materialName = pt.parseString(MATERIAL);
        MaterialPtr material = sceneCache->getMaterial(materialName);

        PrimitivePtr model(new Model(geometry, material));

        if(!model->intersectable()) {
            std::vector<PrimitivePtr> primitives;
            primitives.push_back(model);
            PrimitivePtr aggregate(new BVH(primitives, 1, "equal_count"));
            sceneCache->addPrimitive(name, aggregate);
        } else {
            sceneCache->addPrimitive(name, model);
        }
    }

    static void parseInstance(const PropertyTree& pt, SceneCache* sceneCache) {
        string name = pt.parseString(NAME);
        std::cout << "\ninstance " << name <<std::endl;

        string primitiveName = pt.parseString(MODEL);
        PrimitivePtr primitive = sceneCache->getPrimitive(primitiveName);

        Vector3 position = parseVector3(pt, POSITION);
        std::cout << "-position: " << position << std::endl;
        Quaternion orientation = parseQuaternion(pt, ORIENTATION);
        std::cout << "-orientation: " << orientation << std::endl;
        Vector3 scale = parseVector3(pt, SCALE);
        std::cout << "-scale: " << scale << std::endl;
        Transform toWorld(position, orientation, scale);
        PrimitivePtr instance(new InstancedPrimitive(toWorld, primitive));
        BBox bbox = instance->getAABB();
        std::cout << "BBox min: " << bbox.pMin << std::endl;
        std::cout << "BBox max: " << bbox.pMax << std::endl;

        sceneCache->addInstance(instance);
    }

    static void parseLight(const PropertyTree& pt, SceneCache* sceneCache) {
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
            Vector3 scale = parseVector3(pt, SCALE, Vector3(1.0f, 1.0f, 1.0f));
            string geoName = pt.parseString(GEOMETRY);
            int sampleNum = pt.parseInt(SAMPLE_NUM);
            std::cout << "-radiance: " << radiance << std::endl;
            std::cout << "-position: " << position << std::endl;
            std::cout << "-scale: " << scale << std::endl;
            std::cout << "-orientation: " << orientation << std::endl;
            std::cout << "-geometry: " << geoName << std::endl;
            std::cout << "-samples: " << sampleNum << std::endl;
            GeometryPtr geometry = sceneCache->getGeometry(geoName);
            // TODO: this cause a problem that we can't run time modify
            // the transform for area light since it's not tied in between
            // instance in scene and the transform in area light itself..
            // need to find a way to improve this part
            Transform toWorld(position, orientation, scale);

            AreaLight* areaLight = new AreaLight(radiance, geometry, 
                toWorld, sampleNum);
            light = areaLight;
            // and we need to push this geometry into scene so that it can
            // be intersection tested.....this is gonna be messy for the 
            // current awkward parsing mechanics.......
            TexturePtr white(new ConstantTexture(Color::White));
            MaterialPtr mtl(new LambertMaterial(white));
            PrimitivePtr model(new Model(geometry, mtl, areaLight));
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
            sceneCache->addInstance(instance);
        } else {
            std::cerr << "unrecognized light type " << lightType << std::endl;
            return;
        }
        sceneCache->addLight(light);
    }

    static void parseMaterial(const PropertyTree& pt,
        SceneCache* sceneCache) {
        string materialType = pt.parseString(TYPE);
        string name = pt.parseString(NAME);

        std::cout <<"\n" << materialType <<" material " << name << std::endl;
        MaterialPtr material;
        if(materialType == LAMBERT) {
            string textureName = pt.parseString(DIFFUSE);
            TexturePtr Kd = sceneCache->getTexture(textureName);
            material = MaterialPtr(new LambertMaterial(Kd));
        } else if(materialType == TRANSPARENT) {
            string reflectTextureName = pt.parseString(REFLECTION);
            string refractTextureName = pt.parseString(REFRACTION);
            TexturePtr Kr = sceneCache->getTexture(reflectTextureName);
            TexturePtr Kt = sceneCache->getTexture(refractTextureName);
            float index = pt.parseFloat(REFRACTION_INDEX, 1.5f);
            std::cout << "refraction index: " << index << std::endl;
            material = MaterialPtr(new TransparentMaterial(Kr, Kt, index));
        } else {
            std::cerr << "undefined material type " << 
                materialType << std::endl;
            return;
        }
        sceneCache->addMaterial(name, material);
    }

    static TextureMapping* parseTextureMapping(const PropertyTree& pt) {
        TextureMapping* m;
        string type = pt.parseString(MAPPING);
        std::cout << type << " texture mapping" << std::endl;
        if(type == UV) {
            Vector2 scale = parseVector2(pt, SCALE, Vector2(1.0f, 1.0f));
            Vector2 offset = parseVector2(pt, OFFSET, Vector2::Zero);
            std::cout << "-scale: " << scale << std::endl;
            std::cout << "-offset: " << offset << std::endl;
            m = new UVMapping(scale, offset);
        } else {
            std::cerr << "undefined mapping type " << type << std::endl;
            m = new UVMapping(Vector2(1.0f, 1.0f), Vector2::Zero);
        }
        return m;
    }

    static void parseTexture(const PropertyTree& pt,
        SceneCache* sceneCache) {
        string textureType = pt.parseString(TYPE);
        string name = pt.parseString(NAME);
        std::cout << textureType <<" texture " << name << std::endl;
        TexturePtr texture;
        if(textureType == CONSTANT) {
            Color c = parseColor(pt, COLOR);
            std::cout << "-color " << c << std::endl;
            texture = TexturePtr(new ConstantTexture(c));
        } else if(textureType == IMAGE) {
            TextureMapping* m = parseTextureMapping(pt);
            string filename = 
                sceneCache->resolvePath(pt.parseString(FILENAME));
            float gamma = pt.parseFloat(GAMMA, 1.0f);
            string addressStr = pt.parseString(ADDRESS, REPEAT);
            AddressMode addressMode;
            if(addressStr == REPEAT) {
                addressMode = AddressRepeat;
            } else if(addressStr == CLAMP) {
                addressMode = AddressClamp;
            } else if(addressStr == BORDER) {
                addressMode = AddressBorder;
            }
            std::cout << "-filename: " << filename << std::endl;
            std::cout << "-gamma: " << gamma << std::endl;
            texture = TexturePtr(new ImageTexture(filename, m, 
                addressMode, gamma));
        } else {
            std::cerr << "undefined texture type " << textureType << std::endl;
            return;
        }
        sceneCache->addTexture(name, texture);
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
        path scenePath(filename);

        if(!exists(scenePath) || !pt.read(filename)) {
            return scene;
        }
        SceneCache sceneCache(canonical(scenePath.parent_path()));

        if(setting) {
            parseRenderSetting(pt, setting);
        }
        Film* film = parseFilm(pt, &sceneCache);
        CameraPtr camera = parseCamera(pt, film);

        const PtreeList& nodes = pt.getChildren(); 
        for(PtreeList::const_iterator it = nodes.begin(); 
            it != nodes.end(); it++) {
            std::string key(it->first);
            if(key == MODEL) {
                parseModel(it->second, &sceneCache);
            } else if(key == GEOMETRY) {
                parseGeometry(it->second, &sceneCache);
            } else if(key == INSTANCE) {
                parseInstance(it->second, &sceneCache);
            } else if(key == LIGHT) {
                parseLight(it->second, &sceneCache);
            } else if(key == MATERIAL) {
                parseMaterial(it->second, &sceneCache);
            } else if(key == TEXTURE) {
                parseTexture(it->second, &sceneCache);
            }
        }

        PrimitivePtr aggregate(new BVH(sceneCache.getInstances(),
            1, "equal_count"));
        return ScenePtr(new Scene(aggregate, camera, sceneCache.getLights()));
    }
}
