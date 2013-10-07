#include "GoblinScene.h"

namespace Goblin {

    Scene::Scene(const PrimitivePtr& root, const CameraPtr& camera):
        mAggregate(root), mCamera(camera) {}

    const CameraPtr Scene::getCamera() const {
        return mCamera;
    }

    bool Scene::intersect(const Ray& ray) {
        return mAggregate->intersect(ray);
    }

    void Scene::collectRenderList(RenderList& rList) {
        mAggregate->collectRenderList(rList);
    }
}