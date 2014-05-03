#include "GoblinScene.h"
#include "GoblinTexture.h"
#include "GoblinColor.h"

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
}
