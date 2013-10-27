#include "GoblinScene.h"

namespace Goblin {

    Scene::Scene(const PrimitivePtr& root, const CameraPtr& camera,
        const LightList& lights):
        mAggregate(root), mCamera(camera), mLights(lights) {}

    const CameraPtr Scene::getCamera() const {
        return mCamera;
    }

    const LightList& Scene::getLights() const {
        return mLights;
    }

    bool Scene::intersect(const Ray& ray) {
        return mAggregate->intersect(ray);
    }

    bool Scene::intersect(const Ray& ray, float* epsilon, 
        Intersection* intersection) {
        return mAggregate->intersect(ray, epsilon, intersection);
    }

    void Scene::collectRenderList(RenderList& rList) {
        mAggregate->collectRenderList(rList);
    }
}