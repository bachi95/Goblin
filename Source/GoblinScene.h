#ifndef GOBLIN_SCENE_H
#define GOBLIN_SCENE_H

#include <vector>
#include <boost/shared_ptr.hpp>

namespace Goblin {
    class Camera;
    typedef boost::shared_ptr<Camera> CameraPtr;
    class Model;
    typedef boost::shared_ptr<Model> ModelPtr;
    typedef std::vector<ModelPtr> ModelList;
    class Scene {
    public:
        void addModel(ModelPtr model);
        const ModelList& getModels() const;

        void setCamera(CameraPtr camera);
    private:
        ModelList mModels;
        CameraPtr mCamera;
    };
}

#endif //GOBLIN_SCENE_H