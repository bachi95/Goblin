#ifndef GOBLIN_SCENE_H
#define GOBLIN_SCENE_H

#include "GoblinPrimitive.h"
#include <vector>
#include <boost/shared_ptr.hpp>

namespace Goblin {
    class Camera;
    typedef boost::shared_ptr<Camera> CameraPtr;
    class Ray;

    class Scene {
    public:
        Scene(const PrimitivePtr& root, const CameraPtr& camera);
        const CameraPtr getCamera() const; 
        void collectRenderList(RenderList& rList);
        bool intersect(const Ray& ray);
    private:
        PrimitivePtr mAggregate;
        CameraPtr mCamera;
    };

    typedef boost::shared_ptr<Scene> ScenePtr;
}

#endif //GOBLIN_SCENE_H