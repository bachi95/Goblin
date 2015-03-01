#ifndef GOBLIN_SCENE_H
#define GOBLIN_SCENE_H

#include "GoblinLight.h"
#include "GoblinMaterial.h"
#include "GoblinPrimitive.h"
#include "GoblinTexture.h"
#include "GoblinUtils.h"

#include <boost/filesystem.hpp>

#include <vector>

namespace Goblin {
    class Ray;
    class VolumeRegion;

    class Scene {
    public:
        Scene(const PrimitivePtr& root, const CameraPtr& camera,
            const vector<Light*>& lights, VolumeRegion* volumeRegion);
        ~Scene();
        const CameraPtr getCamera() const; 
        const vector<Light*>& getLights() const;
        const VolumeRegion* getVolumeRegion() const;
        void collectRenderList(RenderList& rList);
        bool intersect(const Ray& ray, IntersectFilter f = NULL) const;
        bool intersect(const Ray& ray, float* epsilon, 
            Intersection* intersection, IntersectFilter f = NULL) const;

        Color evalEnvironmentLight(const Ray& ray) const;
        void getBoundingSphere(Vector3* center, float* radius) const;
    private:
        PrimitivePtr mAggregate;
        CameraPtr mCamera;
        vector<Light*> mLights;
        VolumeRegion* mVolumeRegion;
    };

    using boost::filesystem::path;

    class SceneCache {
    public:
        SceneCache(const path& sceneRoot);
        void addGeometry(const string& name, const Geometry* g);
        void addPrimitive(const string& name, const Primitive* p);
        void addMaterial(const string& name, const MaterialPtr& m);
        void addFloatTexture(const string& name, const FloatTexturePtr& t);
        void addColorTexture(const string& name, const ColorTexturePtr& t);
        void addAreaLight(const string& name, const AreaLight* l);
        void addInstance(const Primitive* i);
        void addLight(Light* l);
        const Geometry* getGeometry(const string& name) const;
        const Primitive* getPrimitive(const string& name) const;
        const MaterialPtr& getMaterial(const string& name) const;
        const FloatTexturePtr& getFloatTexture(const string& name) const;
        const ColorTexturePtr& getColorTexture(const string& name) const;
        const AreaLight* getAreaLight(const string& name) const;
        const PrimitiveList& getInstances() const;
        const vector<Light*>& getLights() const;
        string resolvePath(const string& filename) const;

    private:
        void initDefault();

        typedef map<string, const Geometry*> GeometryMap;
        typedef map<string, const Primitive*> PrimitiveMap;
        typedef map<string, MaterialPtr> MaterialMap;
        typedef map<string, ColorTexturePtr> ColorTextureMap;
        typedef map<string, FloatTexturePtr> FloatTextureMap;
        typedef map<string, const AreaLight*> AreaLightMap;

        GeometryMap mGeometryMap;
        PrimitiveMap mPrimitiveMap;
        MaterialMap mMaterialMap;
        FloatTextureMap mFloatTextureMap;
        ColorTextureMap mColorTextureMap;
        AreaLightMap mAreaLightMap;
        PrimitiveList mInstances;
        vector<Light*> mLights;
        path mSceneRoot;
        string mErrorCode;
    };

}

#endif //GOBLIN_SCENE_H
