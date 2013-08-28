#include "GoblinScene.h"

namespace Goblin {
    void Scene::addModel(ModelPtr model) {
        mModels.push_back(model);
    }

    const ModelList& Scene::getModels() const {
        return mModels;
    }

    void Scene::setCamera(CameraPtr camera) {
        mCamera = camera;
    }
}