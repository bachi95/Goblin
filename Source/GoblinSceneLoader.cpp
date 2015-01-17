#include "GoblinSceneLoader.h"
#include "GoblinBBox.h"
#include "GoblinBVH.h"
#include "GoblinModel.h"
#include "GoblinObjMesh.h"
#include "GoblinPropertyTree.h"
#include "GoblinRenderer.h"
#include "GoblinSphere.h"
#include "GoblinUtils.h"

namespace Goblin {
    static Vector2 parseVector2(const PropertyTree& pt, const char* key,
        const Vector2& fallback = Vector2::Zero) {
        std::vector<float> rv = pt.parseFloatArray(key);
        if(rv.size() != 2) {
            std::cerr << "invalid value for Vector2 " << key << endl;
            return fallback;
        }
        return Vector2(rv[0], rv[1]);
    }

    static Vector3 parseVector3(const PropertyTree& pt, const char* key,
        const Vector3& fallback = Vector3::Zero) {
        std::vector<float> rv = pt.parseFloatArray(key);
        if(rv.size() != 3) {
            std::cerr << "invalid value for Vector3 " << key << endl;
            return fallback;
        }
        return Vector3(rv[0], rv[1], rv[2]);
    }

    static Vector4 parseVector4(const PropertyTree& pt, const char* key,
        const Vector4& fallback = Vector4::Zero) {
        std::vector<float> rv = pt.parseFloatArray(key);
        if(rv.size() != 4) {
            std::cerr << "invalid value for Vector4 " << key << endl;
            return fallback;
        }
        return Vector4(rv[0], rv[1], rv[2], rv[3]);
    }

    static Color parseColor(const PropertyTree& pt, const char* key,
        const Color& fallback = Color::White) {
        std::vector<float> rv = pt.parseFloatArray(key);
        if(rv.size() != 3) {
            std::cerr << "invalid value for Color " << key << endl;
            return fallback;
        }
        return Color(rv[0], rv[1], rv[2], 1.0f);
    }

    static void parseParamSet(const PropertyTree& pt, ParamSet* params) {
        const PtreeList& nodes = pt.getChildren(); 
        cout << "---------------------------------" << endl;
        for(PtreeList::const_iterator it = nodes.begin(); 
            it != nodes.end(); it++) {
            std::string type(it->first);
            const PtreeList& keyValuePairs = it->second.getChildren();
            for(PtreeList::const_iterator kv = keyValuePairs.begin(); 
                kv != keyValuePairs.end(); kv++) {
                std::string key(kv->first);
                cout << type << " " << key << " ";
                if(type == "bool") {
                    bool v = it->second.parseBool(key.c_str());
                    params->setBool(key, v);
                    cout << v << endl; 
                } else if(type == "int") {
                    int v = it->second.parseInt(key.c_str());
                    params->setInt(key, v);
                    cout << v << endl;
                } else if(type == "float") {
                    float v = it->second.parseFloat(key.c_str());
                    params->setFloat(key, v);
                    cout << v << endl;
                } else if(type == "string") {
                    std::string v = it->second.parseString(key.c_str());
                    params->setString(key, v);
                    cout << v << endl;
                }  else if(type == "vec2") {
                    Vector2 v = parseVector2(it->second, key.c_str());
                    params->setVector2(key, v);
                    cout << v << endl;
                }  else if(type == "vec3") {
                    Vector3 v = parseVector3(it->second, key.c_str());
                    params->setVector3(key, v);
                    cout << v << endl;
                }  else if(type == "vec4") {
                    Vector4 v = parseVector4(it->second, key.c_str());
                    params->setVector4(key, v);
                    cout << v << endl;
                } else if(type == "color") {
                    Color v = parseColor(it->second, key.c_str());
                    params->setColor(key, v);
                    cout << v << endl;
                }
            }
        }
        cout << "---------------------------------" << endl;
    }

    static void parseRenderSetting(const PropertyTree& pt, 
        RenderSetting* setting) {
        int samplePerPixel = 1;
        int maxRayDepth = 5;
        int threadNum = boost::thread::hardware_concurrency();
        int bssrdfSampleNum = 4;

        PropertyTree settingPt;
        pt.getChild("render_setting", &settingPt);
        ParamSet settingParams;
        parseParamSet(settingPt, &settingParams);
        samplePerPixel = 
            settingParams.getInt("sample_per_pixel", samplePerPixel);
        maxRayDepth = settingParams.getInt("max_ray_depth", maxRayDepth);
        threadNum = min(settingParams.getInt("thread_num", threadNum), 
            threadNum);
        bssrdfSampleNum = settingParams.getInt("bssrdf_sample_num", 
            bssrdfSampleNum);
        string method = settingParams.getString("render_method", 
            "path_tracing");
        setting->samplePerPixel = samplePerPixel;
        setting->maxRayDepth = maxRayDepth;
        setting->threadNum = threadNum;
        setting->bssrdfSampleNum = bssrdfSampleNum;
        setting->method = method == "path_tracing" ? PathTracing : Whitted;
        cout << "\nrender setting" << endl;
        cout << "-sample per pixel " << samplePerPixel << endl;
        cout << "-thread num " << threadNum << endl;
        cout << "-max ray depth " << maxRayDepth << endl;
        cout << "-bssrdf sample num " << bssrdfSampleNum << endl;
        cout << "-render method " << method << endl;
    }

    SceneLoader::SceneLoader():
        mFilterFactory(new Factory<Filter, const ParamSet&>()),
        mFilmFactory(new Factory<Film, const ParamSet&, Filter*>()),
        mCameraFactory(new Factory<Camera, const ParamSet&, Film*>()),
        mVolumeFactory(new Factory<VolumeRegion, const ParamSet&>()),
        mGeometryFactory(
            new Factory<Geometry, const ParamSet&, const SceneCache&>()),
        mFloatTextureFactory(
            new Factory<Texture<float>, const ParamSet&, const SceneCache&>()),
        mColorTextureFactory(
            new Factory<Texture<Color>, const ParamSet&, const SceneCache&>()),
        mMaterialFactory(
            new Factory<Material, const ParamSet&, const SceneCache&>()),
        mPrimitiveFactory(
            new Factory<Primitive, const ParamSet&, const SceneCache&>()),
        mLightFactory(
            new Factory<Light, const ParamSet&, const SceneCache&>()) {

        // filter
        mFilterFactory->registerCreator("box", new BoxFilterCreator);
        mFilterFactory->registerCreator("triangle", new TriangleFilterCreator);
        mFilterFactory->registerCreator("gaussian", new GaussianFilterCreator);
        mFilterFactory->registerCreator("mitchell", new MitchellFilterCreator);
        mFilterFactory->setDefault("gaussian");
        // film
        mFilmFactory->registerCreator("image", new ImageFilmCreator);
        mFilmFactory->setDefault("image");
        // camera
        mCameraFactory->registerCreator("perspective", 
            new PerspectiveCameraCreator);
        mCameraFactory->setDefault("perspective");
        // volume
        mVolumeFactory->registerCreator("homogeneous", new VolumeCreator);
        mVolumeFactory->setDefault("homogeneous");
        // geometry
        mGeometryFactory->registerCreator("sphere", new SphereGeometryCreator);
        mGeometryFactory->registerCreator("mesh", new MeshGeometryCreator);
        mGeometryFactory->setDefault("sphere");
        // texture
        mFloatTextureFactory->registerCreator("constant", 
            new FloatConstantTextureCreator);
        mFloatTextureFactory->registerCreator("scale", 
            new FloatScaleTextureCreator);
        mFloatTextureFactory->registerCreator("image", 
            new FloatImageTextureCreator);
        mFloatTextureFactory->setDefault("constant");

        mColorTextureFactory->registerCreator("constant", 
            new ColorConstantTextureCreator);
        mColorTextureFactory->registerCreator("scale", 
            new ColorScaleTextureCreator);
        mColorTextureFactory->registerCreator("image", 
            new ColorImageTextureCreator);
        mColorTextureFactory->setDefault("constant");
        // material
        mMaterialFactory->registerCreator("lambert", 
            new LambertMaterialCreator);
        mMaterialFactory->registerCreator("blinn", 
            new BlinnMaterialCreator);
        mMaterialFactory->registerCreator("transparent", 
            new TransparentMaterialCreator);
        mMaterialFactory->registerCreator("mirror", 
            new MirrorMaterialCreator);
        mMaterialFactory->registerCreator("subsurface",
            new SubsurfaceMaterialCreator);
        mMaterialFactory->setDefault("lambert");
        // primitive
        mPrimitiveFactory->registerCreator("model", 
            new ModelPrimitiveCreator);
        mPrimitiveFactory->registerCreator("instance", 
            new InstancePrimitiveCreator);
        mPrimitiveFactory->setDefault("model");
        // light
        mLightFactory->registerCreator("point",
            new PointLightCreator);
        mLightFactory->registerCreator("directional",
            new DirectionalLightCreator);
        mLightFactory->registerCreator("spot",
            new SpotLightCreator);
        mLightFactory->registerCreator("ibl",
            new ImageBasedLightCreator);
        mLightFactory->setDefault("point");
    }

    Filter* SceneLoader::parseFilter(const PropertyTree& pt) {
        cout << "filter" << endl;
        PropertyTree filterPt;
        pt.getChild("filter", &filterPt);
        ParamSet filterParams;
        parseParamSet(filterPt, &filterParams);
        string type = filterParams.getString("type");
        return mFilterFactory->create(type, filterParams);
    }

    Film* SceneLoader::parseFilm(const PropertyTree& pt, Filter* filter) {
        cout << "film" << endl;
        PropertyTree filmPt;
        pt.getChild("film", &filmPt);
        ParamSet filmParams;
        parseParamSet(filmPt, &filmParams);
        string type = filmParams.getString("type");
        return mFilmFactory->create(type, filmParams, filter);
    }

    CameraPtr SceneLoader::parseCamera(const PropertyTree& pt, Film* film) {
        cout << "camera" << endl;
        PropertyTree cameraPt;
        pt.getChild("camera", &cameraPt);
        ParamSet cameraParams;
        parseParamSet(cameraPt, &cameraParams);
        string type = cameraParams.getString("type");
        return CameraPtr(mCameraFactory->create(type, cameraParams, film));
    }

    VolumeRegion* SceneLoader::parseVolume(const PropertyTree& pt) {
        if(!pt.hasChild("volume")) {
            return NULL;
        }
        cout << "volume" << endl;
        PropertyTree volumePt;
        pt.getChild("volume", &volumePt);
        ParamSet volumeParams;
        parseParamSet(volumePt, &volumeParams);
        string type = volumeParams.getString("type");
        return mVolumeFactory->create(type, volumeParams);
    }

    void SceneLoader::parseGeometry(const PropertyTree& pt, 
        SceneCache* sceneCache) {
        cout << "geometry" <<endl;
        ParamSet geometryParams;
        parseParamSet(pt, &geometryParams);
        string type = geometryParams.getString("type");
        string name = geometryParams.getString("name");
        GeometryPtr geometry(mGeometryFactory->create(type, geometryParams,
            *sceneCache));
        geometry->init();
        cout << "vertex num: " << geometry->getVertexNum() << endl;
        cout << "face num: " << geometry->getFaceNum() << endl;
        BBox bbox = geometry->getObjectBound();
        cout << "BBox min: " << bbox.pMin << endl;
        cout << "BBox max: " << bbox.pMax << endl;
        cout << "BBox center: " << 0.5f * (bbox.pMin + bbox.pMax) << endl;

        sceneCache->addGeometry(name, geometry);
    }

    void SceneLoader::parseTexture(const PropertyTree& pt,
        SceneCache* sceneCache) {
        cout << "texture" <<endl;
        ParamSet textureParams;
        parseParamSet(pt, &textureParams);
        string textureFormat = textureParams.getString("format", "color");
        string type = textureParams.getString("type");
        string name = textureParams.getString("name");
        if(textureFormat == "float") {
            FloatTexturePtr texture(mFloatTextureFactory->create(type, 
                textureParams, *sceneCache));
            sceneCache->addFloatTexture(name, texture);
        } else if(textureFormat == "color") {
            ColorTexturePtr texture(mColorTextureFactory->create(type, 
                textureParams, *sceneCache));
            sceneCache->addColorTexture(name, texture);
        } else {
            cerr << "unrecognize texture format" <<
                textureFormat << endl;
        }
    }

    void SceneLoader::parseMaterial(const PropertyTree& pt,
        SceneCache* sceneCache) {
        cout << "material" <<endl;
        ParamSet materialParams;
        parseParamSet(pt, &materialParams);
        string type = materialParams.getString("type");
        string name = materialParams.getString("name");
        MaterialPtr material(mMaterialFactory->create(type, materialParams,
            *sceneCache));
        sceneCache->addMaterial(name, material);
    }

    void SceneLoader::parsePrimitive(const PropertyTree& pt,
        SceneCache* sceneCache) {
        cout << "primitive" <<endl;
        ParamSet primitiveParams;
        parseParamSet(pt, &primitiveParams);
        string type = primitiveParams.getString("type");
        string name = primitiveParams.getString("name");
        PrimitivePtr primitive(mPrimitiveFactory->create(type, primitiveParams,
            *sceneCache));
        sceneCache->addPrimitive(name, primitive);

        if(type == "instance") {
            sceneCache->addInstance(primitive);
        }
    }

    void SceneLoader::parseLight(const PropertyTree& pt, SceneCache* sceneCache,
        int samplePerPixel) {
        cout << "light" << endl;
        ParamSet lightParams;
        parseParamSet(pt, &lightParams);
        lightParams.setInt("sample_per_pixel", samplePerPixel);
        string type = lightParams.getString("type");
        string name = lightParams.getString("name");
        Light* light;
        // TODO figure out a way to get this misfit to generic factory.....
        if(type == "area") {
            Color radiance = lightParams.getColor("radiance");
            Vector3 position = lightParams.getVector3("position");
            Vector4 v = lightParams.getVector4("orientation", 
                Vector4(1.0f, 0.0f, 0.0f, 0.0f));
            Quaternion orientation(v[0], v[1], v[2], v[3]);
            Vector3 scale = lightParams.getVector3("scale", 
                Vector3(1.0f, 1.0f, 1.0f));
            string geoName = lightParams.getString("geometry");
            GeometryPtr geometry = sceneCache->getGeometry(geoName);
            // TODO: this cause a problem that we can't run time modify
            // the transform for area light since it's not tied in between
            // instance in scene and the transform in area light itself..
            // need to find a way to improve this part
            Transform toWorld(position, orientation, scale);
            int sampleNum = lightParams.getInt("sample_num", 1);
            AreaLight* areaLight = new AreaLight(radiance, geometry, 
                toWorld, sampleNum);
            light = areaLight;

            // and we need to push this geometry into scene so that it can
            // be intersection tested.....this is gonna be messy for the 
            // current awkward parsing mechanics.......
            ColorTexturePtr white(new ConstantTexture<Color>(Color::White));
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
            light = mLightFactory->create(type, lightParams,
                *sceneCache);
        }
        sceneCache->addLight(light);
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
        parseRenderSetting(pt, setting);
        Filter* filter = parseFilter(pt);
        Film* film = parseFilm(pt, filter);
        CameraPtr camera = parseCamera(pt, film);
        VolumeRegion* volume = parseVolume(pt);

        PtreeList geometryNodes;
        pt.getChildren("geometry", &geometryNodes);
        for(size_t i = 0; i < geometryNodes.size(); ++i) {
            parseGeometry(geometryNodes[i].second, &sceneCache);
        }
        PtreeList textureNodes;
        pt.getChildren("texture", &textureNodes);
        for(size_t i = 0; i < textureNodes.size(); ++i) {
            parseTexture(textureNodes[i].second, &sceneCache);
        }
        PtreeList materialNodes;
        pt.getChildren("material", &materialNodes);
        for(size_t i = 0; i < materialNodes.size(); ++i) {
            parseMaterial(materialNodes[i].second, &sceneCache);
        }
        PtreeList primitiveNodes;
        pt.getChildren("primitive", &primitiveNodes);
        for(size_t i = 0; i < primitiveNodes.size(); ++i) {
            parsePrimitive(primitiveNodes[i].second, &sceneCache);
        }
        PtreeList lightNodes;
        pt.getChildren("light", &lightNodes);
        for(size_t i = 0; i < lightNodes.size(); ++i) {
            parseLight(lightNodes[i].second, &sceneCache, 
                setting->samplePerPixel);
        }

        PrimitivePtr aggregate(new BVH(sceneCache.getInstances(),
            1, "equal_count"));
        return ScenePtr(new Scene(aggregate, camera, 
            sceneCache.getLights(), volume));
    }
}
