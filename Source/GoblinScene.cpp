#include "GoblinColor.h"
#include "GoblinModel.h"
#include "GoblinParamSet.h"
#include "GoblinSampler.h"
#include "GoblinScene.h"
#include "GoblinSphere.h"
#include "GoblinVolume.h"

namespace Goblin {

    Scene::Scene(const PrimitivePtr& root, const CameraPtr& camera,
        const vector<Light*>& lights, VolumeRegion* volumeRegion):
        mAggregate(root), mCamera(camera), mLights(lights), 
        mVolumeRegion(volumeRegion), mPowerDistribution(NULL) {
        vector<float> lightPowers;
        for (size_t i = 0; i < lights.size(); ++i) {
            lightPowers.push_back(
                lights[i]->power(*this).luminance());
        }
        mPowerDistribution = new CDF1D(lightPowers);
    }

    Scene::~Scene() {        
        for (size_t i = 0; i < mLights.size(); ++i) {
            delete mLights[i];
            mLights[i] = NULL;
        }
        if (mVolumeRegion) {
            delete mVolumeRegion;
            mVolumeRegion = NULL;
        }
        if (mPowerDistribution) {
            delete mPowerDistribution;
            mPowerDistribution = NULL;
        }
        Geometry::clearGeometryCache();
        Primitive::clearAllocatedPrimitives();
        Model::clearRefinedModels();
        ImageTexture<float>::clearImageCache();
        ImageTexture<Color>::clearImageCache();
    }

    const CameraPtr Scene::getCamera() const {
        return mCamera;
    }

    void Scene::getBoundingSphere(Vector3* center, float* radius) const {
        mAggregate->getAABB().getBoundingSphere(center, radius);
    }

    const vector<Light*>& Scene::getLights() const {
        return mLights;
    }

    const VolumeRegion* Scene::getVolumeRegion() const {
        return mVolumeRegion;
    }

    bool Scene::intersect(const Ray& ray, IntersectFilter f) const {
        return mAggregate->intersect(ray, f);
    }

    bool Scene::intersect(const Ray& ray, float* epsilon, 
        Intersection* intersection, IntersectFilter f) const {
        bool isIntersect = mAggregate->intersect(ray, epsilon, intersection, f);
        if (isIntersect) {
            const MaterialPtr& material = intersection->getMaterial();
            material->perturb(&intersection->fragment);
        }
        return isIntersect;
    }

    Color Scene::evalEnvironmentLight(const Ray& ray) const {
        Color Lenv(0.0f);
        for (size_t i = 0; i < mLights.size(); ++i) {
            Lenv += mLights[i]->Le(ray);
        }
        return Lenv;
    }

    void Scene::collectRenderList(RenderList& rList) {
        mAggregate->collectRenderList(rList);
    }

    const Light* Scene::sampleLight(float u, float* pdf) const {
        if (mLights.size() == 0) {
            *pdf = 0.0f;
            return NULL;
        }
        int lightIndex = mPowerDistribution->sampleDiscrete(u, pdf);
        return mLights[lightIndex];
    }

    SceneCache::SceneCache(const path& sceneRoot): 
        mSceneRoot(sceneRoot),
        mErrorCode("error") {
        initDefault();
    }

    void SceneCache::initDefault() {
        Color errorColor = Color::Magenta;
        ColorTexturePtr errorCTexture(new ConstantTexture<Color>(errorColor));
        addColorTexture(mErrorCode, errorCTexture);
        FloatTexturePtr errorFTexture(new ConstantTexture<float>(0.5f));
        addFloatTexture(mErrorCode, errorFTexture); 
        MaterialPtr errorMaterial(new LambertMaterial(errorCTexture));
        addMaterial(mErrorCode, errorMaterial);
        Geometry* errorGeometry = new Sphere(1.0f);
        errorGeometry->init();
        addGeometry(mErrorCode, errorGeometry);
        ParamSet modelParams;
        modelParams.setString("geometry", mErrorCode);
        modelParams.setString("material", mErrorCode);
        const Primitive* errorPrimitive = 
            ModelPrimitiveCreator().create(modelParams, *this);
        addPrimitive(mErrorCode, errorPrimitive);
        addAreaLight(mErrorCode, NULL);
    }

    void SceneCache::addGeometry(const string& name, const Geometry* g) {
        std::pair<string, const Geometry*> pair(name, g);
        mGeometryMap.insert(pair); 
    }

    void SceneCache::addPrimitive(const string& name, const Primitive* p) {
        std::pair<string, const Primitive*> pair(name, p);
        mPrimitiveMap.insert(pair); 
    }

    void SceneCache::addMaterial(const string& name, const MaterialPtr& m) {
        std::pair<string, MaterialPtr> pair(name, m);
        mMaterialMap.insert(pair); 
    }

    void SceneCache::addFloatTexture(const string& name, 
        const FloatTexturePtr& t) {
        std::pair<string, FloatTexturePtr> pair(name, t);
        mFloatTextureMap.insert(pair); 
    }

    void SceneCache::addAreaLight(const string& name, const AreaLight* l) {
        std::pair<string, const AreaLight*> pair(name, l);
        mAreaLightMap.insert(pair); 
    }

    void SceneCache::addColorTexture(const string& name, 
        const ColorTexturePtr& t) {
        std::pair<string, ColorTexturePtr> pair(name, t);
        mColorTextureMap.insert(pair); 
    }

    void SceneCache::addInstance(const Primitive* i) {
        mInstances.push_back(i);
    }

    void SceneCache::addLight(Light* l) {
        mLights.push_back(l);
    }

    const Geometry* SceneCache::getGeometry(const string& name) const {
        GeometryMap::const_iterator it = mGeometryMap.find(name);
        if (it == mGeometryMap.end()) {
            std::cerr << "Geometry " << name << " not defined!\n";
            return mGeometryMap.find(mErrorCode)->second;
        }
        return it->second;
    }

    const Primitive* SceneCache::getPrimitive(const string& name) const {
        PrimitiveMap::const_iterator it = mPrimitiveMap.find(name);
        if (it == mPrimitiveMap.end()) {
            std::cerr << "Primitive " << name << " not defined!\n";
            return mPrimitiveMap.find(mErrorCode)->second;
        }
        return it->second;
    }

    const MaterialPtr& SceneCache::getMaterial(const string& name) const {
        MaterialMap::const_iterator it = mMaterialMap.find(name);
        if (it == mMaterialMap.end()) {
            std::cerr << "Material " << name << " not defined!\n";
            return mMaterialMap.find(mErrorCode)->second;
        }
        return it->second;
    }

    const FloatTexturePtr& SceneCache::getFloatTexture(
        const string& name) const {
        FloatTextureMap::const_iterator it = 
            mFloatTextureMap.find(name);
        if (it == mFloatTextureMap.end()) {
            std::cerr << "Texture " << name << " not defined!\n";
            return mFloatTextureMap.find(mErrorCode)->second;
        }
        return it->second;
    }

    const ColorTexturePtr& SceneCache::getColorTexture(
        const string& name) const {
        ColorTextureMap::const_iterator it = mColorTextureMap.find(name);
        if (it == mColorTextureMap.end()) {
            std::cerr << "Texture " << name << " not defined!\n";
            return mColorTextureMap.find(mErrorCode)->second;
        }
        return it->second;
    }

    const AreaLight* SceneCache::getAreaLight(const string& name) const {
        AreaLightMap::const_iterator it = mAreaLightMap.find(name);
        if (mAreaLightMap.find(name) == mAreaLightMap.end()) {
            std::cerr << "Area Light " << name << " not defined!\n";
            return mAreaLightMap.find(mErrorCode)->second;
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
        if (filePath.is_absolute()) {
            return filePath.generic_string();
        } else {
            return (mSceneRoot / filename).generic_string();
        }
    }
}
