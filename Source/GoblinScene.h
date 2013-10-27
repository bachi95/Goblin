#ifndef GOBLIN_SCENE_H
#define GOBLIN_SCENE_H

#include "GoblinPrimitive.h"
#include "GoblinLight.h"
#include <vector>
#include <boost/shared_ptr.hpp>

namespace Goblin {
    class Camera;
    typedef boost::shared_ptr<Camera> CameraPtr;
    class Ray;

    class Scene {
    public:
        Scene(const PrimitivePtr& root, const CameraPtr& camera,
            const LightList& lights);

        const CameraPtr getCamera() const; 
        const LightList& getLights() const;
        void collectRenderList(RenderList& rList);
        bool intersect(const Ray& ray);
        bool intersect(const Ray& ray, float* epsilon, 
            Intersection* itersection);
    private:
        PrimitivePtr mAggregate;
        LightList mLights;
        CameraPtr mCamera;
    };

    typedef boost::shared_ptr<Scene> ScenePtr;
}

#endif //GOBLIN_SCENE_H