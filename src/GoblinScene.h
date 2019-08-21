#ifndef GOBLIN_SCENE_H
#define GOBLIN_SCENE_H

#include "GoblinBVH.h"
#include "GoblinLight.h"
#include "GoblinMaterial.h"
#include "GoblinPrimitive.h"
#include "GoblinTexture.h"
#include "GoblinUtils.h"

namespace Goblin {
class CDF1D;
class Ray;
class VolumeRegion;

class Scene {
public:
    Scene(const PrimitiveList& inputPrimitives, const CameraPtr& camera,
		std::vector<Geometry*>&& geometries,
		std::vector<Primitive*>&& primitives,
        const std::vector<Light*>& lights, VolumeRegion* volumeRegion);

    ~Scene();

    const CameraPtr getCamera() const;

    const std::vector<Light*>& getLights() const;

    const VolumeRegion* getVolumeRegion() const;

    bool intersect(const Ray& ray, float* epsilon, 
        Intersection* intersection, IntersectFilter f = nullptr) const;

	bool occluded(const Ray& ray, IntersectFilter f = nullptr) const;

    Color evalEnvironmentLight(const Ray& ray) const;

    void getBoundingSphere(Vector3* center, float* radius) const;

    const Light* sampleLight(float u, float* pdf) const;

private:
    BVH mBVH;
    CameraPtr mCamera;
	std::vector<Geometry*> mGeometries;
	std::vector<Primitive*> mPrimitives;
    std::vector<Light*> mLights;
    VolumeRegion* mVolumeRegion;
    CDF1D* mPowerDistribution;
};

class SceneCache {
public:
    SceneCache(const std::string& sceneRoot);
    void addGeometry(const std::string& name, const Geometry* g);
    void addPrimitive(const std::string& name, const Primitive* p);
    void addMaterial(const std::string& name, const MaterialPtr& m);
    void addFloatTexture(const std::string& name, const FloatTexturePtr& t);
    void addColorTexture(const std::string& name, const ColorTexturePtr& t);
    void addAreaLight(const std::string& name, const AreaLight* l);
    void addInstance(const Primitive* i);
    void addLight(Light* l);
    const Geometry* getGeometry(const std::string& name) const;
    const Primitive* getPrimitive(const std::string& name) const;
    const MaterialPtr& getMaterial(const std::string& name) const;
    const FloatTexturePtr& getFloatTexture(const std::string& name) const;
    const ColorTexturePtr& getColorTexture(const std::string& name) const;
    const AreaLight* getAreaLight(const std::string& name) const;
    const PrimitiveList& getInstances() const;
    const std::vector<Light*>& getLights() const;
	std::string resolvePath(const std::string& filename) const;

private:
    void initDefault();

    typedef std::map<std::string, const Geometry*> GeometryMap;
    typedef std::map<std::string, const Primitive*> PrimitiveMap;
    typedef std::map<std::string, MaterialPtr> MaterialMap;
    typedef std::map<std::string, ColorTexturePtr> ColorTextureMap;
    typedef std::map<std::string, FloatTexturePtr> FloatTextureMap;
    typedef std::map<std::string, const AreaLight*> AreaLightMap;

    GeometryMap mGeometryMap;
    PrimitiveMap mPrimitiveMap;
    MaterialMap mMaterialMap;
    FloatTextureMap mFloatTextureMap;
    ColorTextureMap mColorTextureMap;
    AreaLightMap mAreaLightMap;
    PrimitiveList mInstances;
    std::vector<Light*> mLights;
    std::string mSceneRoot;
	std::string mErrorCode;
};

}

#endif //GOBLIN_SCENE_H