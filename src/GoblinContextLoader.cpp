#include "GoblinContextLoader.h"
#include "GoblinBBox.h"
#include "GoblinBDPT.h"
#include "GoblinBVH.h"
#include "GoblinDisk.h"
#include "GoblinLightTracer.h"
#include "GoblinModel.h"
#include "GoblinObjMesh.h"
#include "GoblinAO.h"
#include "GoblinPathtracer.h"
#include "GoblinRenderer.h"
#include "GoblinSphere.h"
#include "GoblinSPPM.h"
#include "GoblinWhitted.h"
#include "GoblinUtils.h"

#include <fstream>
#include "json.hpp"

using json = nlohmann::json;

namespace Goblin {
static const size_t sDelimiterWidth = 75;

static void parseParamSet(const json& jsonContext, ParamSet* params) {
	for (json::const_iterator kv = jsonContext.begin(); kv != jsonContext.end(); kv++) {
		std::string key = kv.key();
		const json& value = kv.value();
		if (value.is_boolean()) {
			params->setBool(key, value);
			std::cout << "bool " << key << " : " << value << std::endl;
		} else if (value.is_number_integer()) {
			params->setInt(key, value);
			std::cout << "int " << key << " : " << value << std::endl;
		} else if (value.is_number_float()) {
			params->setFloat(key, value);
			std::cout << "float " << key << " : " << value << std::endl;
		} else if (value.is_string()) {
			params->setString(key, value);
			std::cout << "string " << key << " : " << value << std::endl;
		} else if (value.is_array()) {
			if (value.size() == 2) {
				Vector2 vec2(value[0], value[1]);
				params->setVector2(key, vec2);
				std::cout << "vec2 " << key << " : " << value << std::endl;
			} else if (value.size() == 3) {
				Vector3 vec3(value[0], value[1], value[2]);
				params->setVector3(key, vec3);
				std::cout << "vec3 " << key << " : " << value << std::endl;
			} else if (value.size() == 4) {
				Vector4 vec4(value[0], value[1], value[2], value[3]);
				params->setVector4(key, vec4);
				std::cout << "vec4 " << key << " : " << value << std::endl;
			}
		}
	}
}

ContextLoader::ContextLoader():
    mFilterFactory(new Factory<Filter, const ParamSet&>()),
    mFilmFactory(new Factory<Film, const ParamSet&, Filter*>()),
    mCameraFactory(new Factory<Camera, const ParamSet&, Film*>()),
    mRendererFactory(new Factory<Renderer, const ParamSet&>()),
    mVolumeFactory(
        new Factory<VolumeRegion, const ParamSet&, const SceneCache&>()),
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
    mCameraFactory->registerCreator("orthographic",
        new OrthographicCameraCreator);
    mCameraFactory->setDefault("perspective");
    // renderer
    mRendererFactory->registerCreator("ao",
        new AORendererCreator);
    mRendererFactory->registerCreator("whitted",
        new WhittedRendererCreator);
    mRendererFactory->registerCreator("path_tracing",
        new PathTracerCreator);
    mRendererFactory->registerCreator("light_tracing",
        new LightTracerCreator);
    mRendererFactory->registerCreator("bdpt",
        new BDPTCreator);
    mRendererFactory->registerCreator("sppm",
        new SPPMCreator);
    mRendererFactory->setDefault("path_tracing");
    // volume
    mVolumeFactory->registerCreator("homogeneous",
        new HomogeneousVolumeCreator);
    mVolumeFactory->registerCreator("heterogeneous",
        new HeterogeneousVolumeCreator);
    mVolumeFactory->setDefault("homogeneous");
    // geometry
    mGeometryFactory->registerCreator("sphere", new SphereGeometryCreator);
    mGeometryFactory->registerCreator("disk", new DiskGeometryCreator);
    mGeometryFactory->registerCreator("mesh", new MeshGeometryCreator);
    mGeometryFactory->setDefault("sphere");
    // texture
    mFloatTextureFactory->registerCreator("constant",
        new FloatConstantTextureCreator);
    mFloatTextureFactory->registerCreator("checkboard",
        new FloatCheckboardTextureCreator);
    mFloatTextureFactory->registerCreator("scale",
        new FloatScaleTextureCreator);
    mFloatTextureFactory->registerCreator("image",
        new FloatImageTextureCreator);
    mFloatTextureFactory->setDefault("constant");

    mColorTextureFactory->registerCreator("constant",
        new ColorConstantTextureCreator);
    mColorTextureFactory->registerCreator("checkboard",
        new ColorCheckboardTextureCreator);
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
    mMaterialFactory->registerCreator("mask",
        new MaskMaterialCreator);
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
    mLightFactory->registerCreator("area",
        new AreaLightCreator);
    mLightFactory->registerCreator("ibl",
        new ImageBasedLightCreator);
    mLightFactory->setDefault("point");
}

static RendererPtr createRenderer(Factory<Renderer, const ParamSet&>& factory,
	const json& jsonContext) {
	std::cout << "render_setting" << std::endl;
	std::cout << std::string(sDelimiterWidth, '-') << std::endl;
	ParamSet setting;
	json::const_iterator it = jsonContext.find("render_setting");
	if (it != jsonContext.end()) {
		parseParamSet(it.value(), &setting);
	}
	std::string method = setting.getString("render_method", "path_tracing");
	std::cout << std::string(sDelimiterWidth, '-') << std::endl;
	return RendererPtr(factory.create(method, setting));
}

static Filter* createFilter(Factory<Filter, const ParamSet&>& factory,
	const json& jsonContext) {
	std::cout << "filter" << std::endl;
	std::cout << std::string(sDelimiterWidth, '-') << std::endl;
	ParamSet filterParams;
	json::const_iterator it = jsonContext.find("filter");
	if (it != jsonContext.end()) {
		parseParamSet(it.value(), &filterParams);
	}
	std::string type = filterParams.getString("type");
	std::cout << std::string(sDelimiterWidth, '-') << std::endl;
	return factory.create(type, filterParams);
}

static Film* createFilm(Factory<Film, const ParamSet&, Filter*>& factory,
	const json& jsonContext, Filter* filter) {
	std::cout << "film" << std::endl;
	std::cout << std::string(sDelimiterWidth, '-') << std::endl;
	ParamSet filmParams;
	json::const_iterator it = jsonContext.find("film");
	if (it != jsonContext.end()) {
		parseParamSet(it.value(), &filmParams);
	}
	std::string type = filmParams.getString("type");
	std::cout << std::string(sDelimiterWidth, '-') << std::endl;
	return factory.create(type, filmParams, filter);
}

static CameraPtr createCamera(Factory<Camera, const ParamSet&, Film*>& cameraFactory,
	Factory<Primitive, const ParamSet&, const SceneCache&>& primitiveFactory,
	const json& jsonContext, Film* film, SceneCache* sceneCache) {
	std::cout << "camera" << std::endl;
	std::cout << std::string(sDelimiterWidth, '-') << std::endl;
	ParamSet cameraParams;
	json::const_iterator it = jsonContext.find("camera");
	if (it != jsonContext.end()) {
		parseParamSet(it.value(), &cameraParams);
	}
	std::string type = cameraParams.getString("type");
	// need to add lens into scene so that it can be intersected by
	// light particles
	float lensRadius = cameraParams.getFloat("lens_radius");
	if (lensRadius != 0.0f) {
		// we need to push this geometry into scene so that it can
		// be intersection tested. the run time material/model
		// creation is pretty awkward at this moment...definitly should
		// be improved....
		ParamSet modelParams;
		Geometry* lensGeom = new Disk(lensRadius);
		std::string lensGeomName = type + "_lens_geom";
		sceneCache->addGeometry(lensGeomName, lensGeom);
		modelParams.setString("geometry", lensGeomName);
		modelParams.setBool("is_camera_lens", true);
		ColorTexturePtr black(new ConstantTexture<Color>(Color::Black));
		MaterialPtr mtl(new LambertMaterial(black));
		std::string materialName = type + "_lens_material";
		sceneCache->addMaterial(materialName, mtl);
		modelParams.setString("material", materialName);
		const Primitive* model(primitiveFactory.create("model",
			modelParams, *sceneCache));
		std::string modelName = type + "_lens_model";
		sceneCache->addPrimitive(modelName, model);
		cameraParams.setString("model", modelName);
		const Primitive* instance(primitiveFactory.create("instance",
			cameraParams, *sceneCache));
		sceneCache->addInstance(instance);
	}
	std::cout << std::string(sDelimiterWidth, '-') << std::endl;
	return CameraPtr(cameraFactory.create(type, cameraParams, film));
}

static VolumeRegion* createVolume(
	Factory<VolumeRegion, const ParamSet&, const SceneCache&>& factory,
	const json& jsonContext, SceneCache* sceneCache) {
	std::cout << "volume" << std::endl;
	std::cout << std::string(sDelimiterWidth, '-') << std::endl;
	ParamSet volumeParams;
	json::const_iterator it = jsonContext.find("volume");
	if (it != jsonContext.end()) {
		parseParamSet(it.value(), &volumeParams);
	} else {
		return nullptr;
	}
	std::string type = volumeParams.getString("type");
	return factory.create(type, volumeParams, *sceneCache);
}

static void createGeometries(
	Factory<Geometry, const ParamSet&, const SceneCache&>& factory,
	const json& jsonContext, SceneCache* sceneCache) {
	json::const_iterator it = jsonContext.find("geometries");
	if (it != jsonContext.end()) {
		const json& geometriesList = it.value();
		assert(geometriesList.is_array());
		for (size_t i = 0; i < geometriesList.size(); ++i) {
			const json& geometryContext = geometriesList[i];
			assert(geometryContext.is_object());
			std::cout << "geometry" << std::endl;
			std::cout << std::string(sDelimiterWidth, '-') << std::endl;
			ParamSet geometryParams;
			parseParamSet(geometryContext, &geometryParams);
			std::string type = geometryParams.getString("type");
			std::string name = geometryParams.getString("name");
			Geometry* geometry(factory.create(type, geometryParams,
				*sceneCache));
			geometry->init();
			std::cout << "vertex num: " << geometry->getVertexNum() << std::endl;
			std::cout << "face num: " << geometry->getFaceNum() << std::endl;
			BBox bbox = geometry->getObjectBound();
			std::cout << "BBox min: " << bbox.pMin << std::endl;
			std::cout << "BBox max: " << bbox.pMax << std::endl;
			std::cout << "BBox center: " << bbox.center() << std::endl;
			std::cout << std::string(sDelimiterWidth, '-') << std::endl;
			sceneCache->addGeometry(name, geometry);
		}
	}
}

static void createTextures(
	Factory<Texture<float>, const ParamSet&, const SceneCache&>& floatTextureFactory,
	Factory<Texture<Color>, const ParamSet&, const SceneCache&>& colorTextureFactory,
	const json& jsonContext, SceneCache* sceneCache) {
	json::const_iterator it = jsonContext.find("textures");
	if (it != jsonContext.end()) {
		const json& texturesList = it.value();
		assert(texturesList.is_array());
		for (size_t i = 0; i < texturesList.size(); ++i) {
			const json& textureContext = texturesList[i];
			assert(textureContext.is_object());
			std::cout << "texture" << std::endl;
			std::cout << std::string(sDelimiterWidth, '-') << std::endl;
			ParamSet textureParams;
			parseParamSet(textureContext, &textureParams);
			std::string type = textureParams.getString("type");
			std::string name = textureParams.getString("name");
			std::string textureFormat = textureParams.getString("format", "color");
			if (textureFormat == "float") {
				FloatTexturePtr texture(floatTextureFactory.create(type,
					textureParams, *sceneCache));
				sceneCache->addFloatTexture(name, texture);
			} else if (textureFormat == "color") {
				ColorTexturePtr texture(colorTextureFactory.create(type,
					textureParams, *sceneCache));
				sceneCache->addColorTexture(name, texture);
			} else {
				std::cerr << "unrecognize texture format" <<
					textureFormat << std::endl;
			}
		}
	}
}

static void createMaterials(
	Factory<Material, const ParamSet&, const SceneCache&>& factory,
	const json& jsonContext, SceneCache* sceneCache) {
	json::const_iterator it = jsonContext.find("materials");
	if (it != jsonContext.end()) {
		const json& materialsList = it.value();
		assert(materialsList.is_array());
		for (size_t i = 0; i < materialsList.size(); ++i) {
			const json& materialContext = materialsList[i];
			assert(materialContext.is_object());
			std::cout << "material" << std::endl;
			std::cout << std::string(sDelimiterWidth, '-') << std::endl;
			ParamSet materialParams;
			parseParamSet(materialContext, &materialParams);
			std::string type = materialParams.getString("type");
			std::string name = materialParams.getString("name");
			MaterialPtr material(factory.create(type, materialParams,
				*sceneCache));
			sceneCache->addMaterial(name, material);
			std::cout << std::string(sDelimiterWidth, '-') << std::endl;
		}
	}
}

static void createPrimitives(
	Factory<Primitive, const ParamSet&, const SceneCache&>& factory,
	const json& jsonContext, SceneCache* sceneCache) {
	json::const_iterator it = jsonContext.find("primitives");
	if (it != jsonContext.end()) {
		const json& primitivesList = it.value();
		assert(primitivesList.is_array());
		for (size_t i = 0; i < primitivesList.size(); ++i) {
			const json& primitiveContext = primitivesList[i];
			assert(primitiveContext.is_object());
			std::cout << "primitive" << std::endl;
			std::cout << std::string(sDelimiterWidth, '-') << std::endl;
			ParamSet primitiveParams;
			parseParamSet(primitiveContext, &primitiveParams);
			std::string type = primitiveParams.getString("type");
			std::string name = primitiveParams.getString("name");
			const Primitive* primitive = factory.create(type,
				primitiveParams, *sceneCache);
			sceneCache->addPrimitive(name, primitive);
			BBox bbox = primitive->getAABB();
			std::cout << "BBox min: " << bbox.pMin << std::endl;
			std::cout << "BBox max: " << bbox.pMax << std::endl;
			std::cout << "BBox center: " << bbox.center() << std::endl;
			if (type == "instance") {
				sceneCache->addInstance(primitive);
			}
			std::cout << std::string(sDelimiterWidth, '-') << std::endl;
		}
	}
}

static void createLights(
	Factory<Light, const ParamSet&, const SceneCache&>& lightFactory,
	Factory<Primitive, const ParamSet&, const SceneCache&>& primitiveFactory,
	const json& jsonContext, SceneCache* sceneCache) {
	json::const_iterator it = jsonContext.find("lights");
	if (it != jsonContext.end()) {
		const json& lightsList = it.value();
		assert(lightsList.is_array());
		for (size_t i = 0; i < lightsList.size(); ++i) {
			const json& lightContext = lightsList[i];
			assert(lightContext.is_object());
			std::cout << "light" << std::endl;
			std::cout << std::string(sDelimiterWidth, '-') << std::endl;
			ParamSet lightParams;
			parseParamSet(lightContext, &lightParams);
			std::string type = lightParams.getString("type");
			std::string name = lightParams.getString("name");

			Light* light = lightFactory.create(type, lightParams, *sceneCache);
			sceneCache->addLight(light);
			if (type == "area") {
				// we need to push this geometry into scene so that it can
				// be intersection tested. the run time material/model
				// creation is pretty awkward at this moment...definitly should
				// be improved....
				ParamSet modelParams;
				modelParams.setString("geometry",
					lightParams.getString("geometry"));
				ColorTexturePtr black(new ConstantTexture<Color>(Color::Black));
				MaterialPtr mtl(new LambertMaterial(black));
				string materialName = type + "_" + name + "_material";
				sceneCache->addMaterial(materialName, mtl);
				modelParams.setString("material", materialName);
				const AreaLight* areaLight = static_cast<AreaLight*>(light);
				sceneCache->addAreaLight(name, areaLight);
				modelParams.setString("area_light", name);
				const Primitive* model(primitiveFactory.create("model",
					modelParams, *sceneCache));
				string modelName = type + "_" + name + "_model";
				sceneCache->addPrimitive(modelName, model);
				lightParams.setString("model", modelName);
				const Primitive* instance(primitiveFactory.create("instance",
					lightParams, *sceneCache));
				sceneCache->addInstance(instance);
			}
			std::cout << std::string(sDelimiterWidth, '-') << std::endl;
		}
	}
}

RenderContext* ContextLoader::load(const string& filename) {
	std::ifstream jsonFileStream(filename);
	if (!jsonFileStream.is_open()) {
		std::cerr << "error reading scene file: " << filename << std::endl;
		return nullptr;
	}
	json jsonContext;
	try {
		jsonContext = json::parse(jsonFileStream);
	} catch (json::parse_error& e) {
		std::cout << e.what() << std::endl;
		return nullptr;
	}

	// TODO replace this workaround as soon as we have std::filesystem
	std::string sceneDir = ".";
	std::size_t substrIndex = filename.find_last_of('/');
	if (substrIndex != std::string::npos) {
		sceneDir = filename.substr(0, substrIndex);
	} else {
		substrIndex = filename.find_last_of('\\');
		if (substrIndex != std::string::npos) {
			sceneDir = filename.substr(0, substrIndex);
		}
	}
    SceneCache sceneCache(sceneDir);

	RendererPtr renderer = createRenderer(*mRendererFactory, jsonContext);
	Filter* filter = createFilter(*mFilterFactory, jsonContext);
	Film* film = createFilm(*mFilmFactory, jsonContext, filter);
	CameraPtr camera = createCamera(*mCameraFactory, *mPrimitiveFactory,
		jsonContext, film, &sceneCache);
	VolumeRegion* volume = createVolume(*mVolumeFactory,
		jsonContext, &sceneCache);
	createGeometries(*mGeometryFactory, jsonContext, &sceneCache);
	createTextures(*mFloatTextureFactory, *mColorTextureFactory,
		jsonContext, &sceneCache);
	createMaterials(*mMaterialFactory, jsonContext, &sceneCache);
	createPrimitives(*mPrimitiveFactory, jsonContext, &sceneCache);
	createLights(*mLightFactory, *mPrimitiveFactory, jsonContext, &sceneCache);

    PrimitivePtr aggregate(new BVH(sceneCache.getInstances(),
        1, "equal_count"));
    ScenePtr scene(new Scene(aggregate, camera,
        sceneCache.getLights(), volume));

    RenderContext* ctx = new RenderContext(renderer, scene);
    return ctx;
}

} // namespace Goblin