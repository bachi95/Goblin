#ifndef GOBLIN_CAMERA_H
#define GOBLIN_CAMERA_H

#include "GoblinMatrix.h"
#include "GoblinVector.h"

namespace Goblin {
    class Camera {
    public:
        Camera();
        ~Camera();

        Vector3& position();
        Matrix4 view() const;
        Matrix4 proj() const;

        void setLens(float fovY, float aspect, float nearm, float far);
        void strafe(float d);
        void walk(float d);
        void pitch(float angle);
        void rotateY(float angle);

        void rebuildView();
    private:
        Vector3 mPosition;
        Vector3 mRight;
        Vector3 mUp;
        Vector3 mLook;

        Matrix4 mView;
        Matrix4 mProj;
    };
}

#endif //GOBLIN_CAMERA_H