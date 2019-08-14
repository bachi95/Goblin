#include "GoblinContextLoader.h"
#include "GoblinBBox.h"
#include "GoblinBDPT.h"
#include "GoblinBVH.h"
#include "GoblinCamera.h"
#include "GoblinDisk.h"
#include "GoblinFilter.h"
#include "GoblinFilm.h"
#include "GoblinLightTracer.h"
#include "GoblinModel.h"
#include "GoblinPolygonMesh.h"
#include "GoblinAO.h"
#include "GoblinPathtracer.h"
#include "GoblinPrimitive.h"
#include "GoblinRenderer.h"
#include "GoblinRenderContext.h"
#include "GoblinScene.h"
#include "GoblinSphere.h"
#include "GoblinSPPM.h"
#include "GoblinTexture.h"
#include "GoblinUtils.h"
#include "GoblinVolume.h"
#include "GoblinWhitted.h"

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

static RendererPtr createRenderer(const json& jsonContext) {
	std::cout << "render_setting" << std::endl;
	std::cout << std::string(sDelimiterWidth, '-') << std::endl;
	ParamSet setting;
	json::const_iterator it = jsonContext.find("render_setting");
	if (it != jsonContext.end()) {
		parseParamSet(it.value(), &setting);
	}
	std::string method = setting.getString("render_method", "path_tracing");
	std::cout << std::string(sDelimiterWidth, '-') << std::endl;
	if (method == "ao") {
		return RendererPtr(createAO(setting));
	} else if (method == "whitted") {
		return RendererPtr(createWhitted(setting));
	} else if (method == "path_tracing") {
		return RendererPtr(createPathTracer(setting));
	} else if (method == "light_tracing") {
		return RendererPtr(createLightTracer(setting));
	} else if (method == "bdpt") {
		return RendererPtr(createBDPT(setting));
	} else if (method == "sppm") {
		return RendererPtr(createSPPM(setting));
	} else {
		return RendererPtr(createPathTracer(setting));
	}
}

static Filter* createFilter(const json& jsonContext) {
	std::cout << "filter" << std::endl;
	std::cout << std::string(sDelimiterWidth, '-') << std::endl;
	ParamSet filterParams;
	json::const_iterator it = jsonContext.find("filter");
	if (it != jsonContext.end()) {
		parseParamSet(it.value(), &filterParams);
	}
	std::string type = filterParams.getString("type");
	std::cout << std::string(sDelimiterWidth, '-') << std::endl;
	if (type == "box") {
		return createBoxFilter(filterParams);
	} else if (type == "triangle") {
		return createTriangleFilter(filterParams);
	} else if (type == "gaussian") {
		return createGaussianFilter(filterParams);
	} else if (type == "mitchell") {
		return createMitchellFilter(filterParams);
	} else {
		return createGaussianFilter(filterParams);
	}
}

static Film* createFilm(const json& jsonContext,
	const std::string& defaultOutputPath) {
	std::cout << "film" << std::endl;
	std::cout << std::string(sDelimiterWidth, '-') << std::endl;
	ParamSet filmParams;
	json::const_iterator it = jsonContext.find("film");
	if (it != jsonContext.end()) {
		parseParamSet(it.value(), &filmParams);
	}
	if (!filmParams.hasString("file")) {
		filmParams.setString("file", defaultOutputPath);
	}
	std::cout << std::string(sDelimiterWidth, '-') << std::endl;
	return createImageFilm(filmParams, createFilter(jsonContext));
}

static CameraPtr createCamera(const json& jsonContext, SceneCache* sceneCache,
	const std::string& defaultOutputPath) {
	std::cout << "camera" << std::endl;
	std::cout << std::string(sDelimiterWidth, '-') << std::endl;
	ParamSet cameraParams;
	json cameraContext;
	json::const_iterator it = jsonContext.find("camera");
	if (it != jsonContext.end()) {
		cameraContext = it.value();
		parseParamSet(cameraContext, &cameraParams);
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
		const Primitive* model(createModel(modelParams, *sceneCache));
		std::string modelName = type + "_lens_model";
		sceneCache->addPrimitive(modelName, model);
		cameraParams.setString("model", modelName);
		const Primitive* instance(createInstance(cameraParams, *sceneCache));
		sceneCache->addInstance(instance);
	}
	std::cout << std::string(sDelimiterWidth, '-') << std::endl;
	if (type == "perspective") {
		return CameraPtr(createPerspectiveCamera(cameraParams,
			createFilm(cameraContext, defaultOutputPath)));
	} else if (type == "orthographic") {
		return CameraPtr(createOrthographicCamera(cameraParams,
			createFilm(cameraContext, defaultOutputPath)));
	} else {
		return CameraPtr(createPerspectiveCamera(cameraParams,
			createFilm(cameraContext, defaultOutputPath)));
	}
}

static VolumeRegion* createVolume(const json& jsonContext,
	SceneCache* sceneCache) {
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
	if (type == "homogeneous") {
		return createHomogeneousVolume(volumeParams, *sceneCache);
	} else if (type == "heterogeneous") {
		return createHeterogeneousVolume(volumeParams, *sceneCache);
	} else {
		return createHomogeneousVolume(volumeParams, *sceneCache);
	}
}

static void createGeometries(const json& jsonContext, SceneCache* sceneCache) {
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
			Geometry* geometry = nullptr;
			if (type == "sphere") {
				geometry = createSphere(geometryParams, *sceneCache);
			} else if (type == "mesh") {
				geometry = createPolygonMesh(geometryParams, *sceneCache);
			} else if (type == "disk") {
				geometry = createDisk(geometryParams, *sceneCache);
			} else {
				geometry = createSphere(geometryParams, *sceneCache);
			}
			BBox bbox = geometry->getObjectBound();
			std::cout << "BBox min: " << bbox.pMin << std::endl;
			std::cout << "BBox max: " << bbox.pMax << std::endl;
			std::cout << "BBox center: " << bbox.center() << std::endl;
			std::cout << std::string(sDelimiterWidth, '-') << std::endl;
			sceneCache->addGeometry(name, geometry);
		}
	}
}

static void createTextures(const json& jsonContext, SceneCache* sceneCache) {
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
				FloatTexturePtr texture;
				if (type == "constant") {
					texture.reset(createFloatConstantTexture(textureParams));
				} else if (type == "checkerboard") {
					texture.reset(createFloatCheckerboardTexture(
						textureParams, *sceneCache));
				} else if (type == "scale") {
					texture.reset(createFloatScaleTexture(
						textureParams, *sceneCache));
				} else if (type == "image") {
					texture.reset(createFloatImageTexture(
						textureParams, *sceneCache));
				} else {
					texture.reset(createFloatConstantTexture(
						textureParams));
				}
				sceneCache->addFloatTexture(name, texture);
			} else if (textureFormat == "color") {
				ColorTexturePtr texture;
				if (type == "constant") {
					texture.reset(createColorConstantTexture(textureParams));
				}
				else if (type == "checkerboard") {
					texture.reset(createColorCheckerboardTexture(
						textureParams, *sceneCache));
				}
				else if (type == "scale") {
					texture.reset(createColorScaleTexture(
						textureParams, *sceneCache));
				}
				else if (type == "image") {
					texture.reset(createColorImageTexture(
						textureParams, *sceneCache));
				}
				else {
					texture.reset(createColorConstantTexture(
						textureParams));
				}
				sceneCache->addColorTexture(name, texture);
			} else {
				std::cerr << "unrecognize texture format" <<
					textureFormat << std::endl;
			}
		}
	}
}

static void createMaterials(const json& jsonContext, SceneCache* sceneCache) {
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
			MaterialPtr material;
			if (type == "lambert") {
				material.reset(createLambertMaterial(materialParams,
					*sceneCache));
			} else if (type == "blinn") {
				material.reset(createBlinnMaterial(materialParams,
					*sceneCache));
			} else if (type == "transparent") {
				material.reset(createTransparentMaterial(materialParams,
					*sceneCache));
			} else if (type == "mirror") {
				material.reset(createMirrorMaterial(materialParams,
					*sceneCache));
			} else if (type == "subsurface") {
				material.reset(createSubsurfaceMaterial(materialParams,
					*sceneCache));
			} else if (type == "mask") {
				material.reset(createMaskMaterial(materialParams,
					*sceneCache));
			} else {
				material.reset(createLambertMaterial(materialParams,
					*sceneCache));
			}
			sceneCache->addMaterial(name, material);
			std::cout << std::string(sDelimiterWidth, '-') << std::endl;
		}
	}
}

static void createPrimitives(const json& jsonContext, SceneCache* sceneCache) {
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
			Primitive* primitive = nullptr;
			if (type == "model") {
				primitive = createModel(primitiveParams, *sceneCache);
			} else if (type == "instance") {
				primitive = createInstance(primitiveParams, *sceneCache);
			} else {
				primitive = createModel(primitiveParams, *sceneCache);
			}
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

static void createLights(const json& jsonContext, SceneCache* sceneCache) {
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

			Light* light = nullptr;
			if (type == "point") {
				light = createPointLight(lightParams, *sceneCache);
			} else if (type == "directional") {
				light = createDirectionalLight(lightParams, *sceneCache);
			} else if (type == "spot") {
				light = createSpotLight(lightParams, *sceneCache);
			} else if (type == "area") {
				light = createAreaLight(lightParams, *sceneCache);
			} else if (type == "ibl") {
				light = createImageBasedLight(lightParams, *sceneCache);
			} else {
				light = createPointLight(lightParams, *sceneCache);
			}
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
				std::string materialName = type + "_" + name + "_material";
				sceneCache->addMaterial(materialName, mtl);
				modelParams.setString("material", materialName);
				const AreaLight* areaLight = static_cast<AreaLight*>(light);
				sceneCache->addAreaLight(name, areaLight);
				modelParams.setString("area_light", name);
				const Primitive* model(createModel(modelParams, *sceneCache));
				std::string modelName = type + "_" + name + "_model";
				sceneCache->addPrimitive(modelName, model);
				lightParams.setString("model", modelName);
				const Primitive* instance(createInstance(lightParams, *sceneCache));
				sceneCache->addInstance(instance);
			}
			std::cout << std::string(sDelimiterWidth, '-') << std::endl;
		}
	}
}

RenderContext* ContextLoader::load(const std::string& filename) {
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
	std::size_t extensionIndex = filename.find_last_of('.');
	// figure out a default output path
	std::string defaultOutputPath;
	if (extensionIndex != std::string::npos &&
		extensionIndex < (filename.length() - 1) &&
		filename[extensionIndex + 1] != '/' &&
		filename[extensionIndex + 1] != '\\') {
		defaultOutputPath = filename.substr(0, extensionIndex) +
			std::string(".exr");
	} else {
		defaultOutputPath = filename + std::string(".exr");
	}

	RendererPtr renderer = createRenderer(jsonContext);
	CameraPtr camera = createCamera(
		jsonContext, &sceneCache, defaultOutputPath);
	VolumeRegion* volume = createVolume(jsonContext, &sceneCache);
	createGeometries(jsonContext, &sceneCache);
	createTextures(jsonContext, &sceneCache);
	createMaterials(jsonContext, &sceneCache);
	createPrimitives(jsonContext, &sceneCache);
	createLights(jsonContext, &sceneCache);

    ScenePtr scene(new Scene(sceneCache.getInstances(), camera,
        sceneCache.getLights(), volume));

    RenderContext* ctx = new RenderContext(renderer, scene);
    return ctx;
}

} // namespace Goblin