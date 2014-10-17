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
        bool intersect(const Ray& ray);
        bool intersect(const Ray& ray, float* epsilon, 
            Intersection* itersection);
        void getBoundingSphere(Vector3* center, float* radius) const;
    private:
        PrimitivePtr mAggregate;
        CameraPtr mCamera;
        vector<Light*> mLights;
        VolumeRegion* mVolumeRegion;
    };

    using boost::filesystem::path;
    struct RenderSetting;

    class SceneCache {
    public:
        SceneCache(const path& sceneRoot);
        void addGeometry(const string& name, const GeometryPtr& g);
        void addPrimitive(const string& name, const PrimitivePtr& p);
        void addMaterial(const string& name, const MaterialPtr& m);
        void addFloatTexture(const string& name, const FloatTexturePtr& t);
        void addColorTexture(const string& name, const ColorTexturePtr& t);
        void addInstance(const PrimitivePtr& i);
        void addLight(Light* l);
        const GeometryPtr& getGeometry(const string& name) const;
        const PrimitivePtr& getPrimitive(const string& name) const;
        const MaterialPtr& getMaterial(const string& name) const;
        const FloatTexturePtr& getFloatTexture(const string& name) const;
        const ColorTexturePtr& getColorTexture(const string& name) const;
        const PrimitiveList& getInstances() const;
        const vector<Light*>& getLights() const;
        string resolvePath(const string& filename) const;

    private:
        void initDefault();

        typedef map<string, GeometryPtr> GeometryMap;
        typedef map<string, PrimitivePtr> PrimitiveMap;
        typedef map<string, MaterialPtr> MaterialMap;
        typedef map<string, ColorTexturePtr> ColorTextureMap;
        typedef map<string, FloatTexturePtr> FloatTextureMap;

        GeometryMap mGeometryMap;
        PrimitiveMap mPrimitiveMap;
        MaterialMap mMaterialMap;
        FloatTextureMap mFloatTextureMap;
        ColorTextureMap mColorTextureMap;
        PrimitiveList mInstances;
        vector<Light*> mLights;
        path mSceneRoot;
        string mErrorCode;
    };

}

#endif //GOBLIN_SCENE_H
