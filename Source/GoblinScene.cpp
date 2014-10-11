#include "GoblinColor.h"
#include "GoblinModel.h"
#include "GoblinScene.h"
#include "GoblinSphere.h"

namespace Goblin {

    Scene::Scene(const PrimitivePtr& root, const CameraPtr& camera,
        const vector<Light*>& lights):
        mAggregate(root), mCamera(camera), mLights(lights) {}

    Scene::~Scene() {        
        for(size_t i = 0; i < mLights.size(); ++i) {
            delete mLights[i];
            mLights[i] = NULL;
        }
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

    bool Scene::intersect(const Ray& ray) {
        return mAggregate->intersect(ray);
    }

    bool Scene::intersect(const Ray& ray, float* epsilon, 
        Intersection* intersection) {
        bool isIntersect = mAggregate->intersect(ray, epsilon, intersection);
        if(isIntersect) {
            const MaterialPtr& material = intersection->getMaterial();
            material->perturb(&intersection->fragment);
        }
        return isIntersect;
    }

    void Scene::collectRenderList(RenderList& rList) {
        mAggregate->collectRenderList(rList);
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
        GeometryPtr errorGeometry(new Sphere(1.0f));
        addGeometry(mErrorCode, errorGeometry);
        PrimitivePtr errorPrimitive(new Model(errorGeometry, errorMaterial));
        addPrimitive(mErrorCode, errorPrimitive);
    }

    void SceneCache::addGeometry(const string& name, const GeometryPtr& g) {
        std::pair<string, GeometryPtr> pair(name, g);
        mGeometryMap.insert(pair); 
    }

    void SceneCache::addPrimitive(const string& name, const PrimitivePtr& p) {
        std::pair<string, PrimitivePtr> pair(name, p);
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

    void SceneCache::addColorTexture(const string& name, 
        const ColorTexturePtr& t) {
        std::pair<string, ColorTexturePtr> pair(name, t);
        mColorTextureMap.insert(pair); 
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

    const FloatTexturePtr& SceneCache::getFloatTexture(
        const string& name) const {
        FloatTextureMap::const_iterator it = 
            mFloatTextureMap.find(name);
        if(it == mFloatTextureMap.end()) {
            std::cerr << "Texture " << name << " not defined!\n";
            return mFloatTextureMap.find(mErrorCode)->second;
        }
        return it->second;
    }

    const ColorTexturePtr& SceneCache::getColorTexture(
        const string& name) const {
        ColorTextureMap::const_iterator it = mColorTextureMap.find(name);
        if(it == mColorTextureMap.end()) {
            std::cerr << "Texture " << name << " not defined!\n";
            return mColorTextureMap.find(mErrorCode)->second;
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
}
