#ifndef GOBLIN_PATH_VERTEX_H
#define GOBLIN_PATH_VERTEX_H

#include "GoblinColor.h"
#include "GoblinCamera.h"
#include "GoblinPrimitive.h"
#include "GoblinLight.h"
#include "GoblinUtils.h"

namespace Goblin {
    class PathVertex{
    public:
        PathVertex(): throughput(0.0f),
            light(NULL), material(NULL), pdfForward(0.0f), pdfBackward(0.0f),
            isCameraLens(false), isSpecular(false), G(0.0f) {}

        PathVertex(const Color& t, const Vector3& p, const Vector3& n,
            const Light* l, float pForward = 0.0f, float pBackward = 0.0f):
            throughput(t),
            fragment(p, n, Vector2::Zero, Vector3::Zero, Vector3::Zero),
            light(l), material(NULL),
            pdfForward(pForward), pdfBackward(pBackward),
            isCameraLens(false), isSpecular(false), G(0.0f) {}

        PathVertex(const Color& t, const Vector3& p, const Vector3& n,
            const Camera* c, float pForward = 0.0f,
            float pBackward = 0.0f):
            throughput(t),
            fragment(p, n, Vector2::Zero, Vector3::Zero, Vector3::Zero),
            light(NULL), material(NULL),
            pdfForward(pForward), pdfBackward(pBackward),
            isCameraLens(true), isSpecular(false), G(0.0f) {}

        PathVertex(const Color& t, const Intersection& isect,
            float pForward = 0.0f, float pBackward = 0.0f, bool spec = false):
            throughput(t), fragment(isect.fragment), light(isect.getLight()),
            material(isect.getMaterial().get()),
            pdfForward(pForward), pdfBackward(pBackward),
            isCameraLens(isect.isCameraLens()), isSpecular(spec), G(0.0f) {}

        const Vector3& getPosition() const {
            return fragment.getPosition();
        }

        const Vector3& getNormal() const {
            return fragment.getNormal();
        }

        bool isLight() const {
            return light != NULL;
        }

        const Light* getLight() const {
            return light;
        }

        Color throughput;
        Fragment fragment;
        const Light* light;
        const Material* material;
        float pdfForward;
        float pdfBackward;
        bool isCameraLens;
        bool isSpecular;
        float G;
    };

    class BDPTMISNode {
    public:
        BDPTMISNode(float towardLight = 0.0f, float towardEye = 0.0f,
            bool spec = false): pTowardLight(towardLight),
            pTowardEye(towardEye), isSpecular(spec) {}
        float pTowardLight;
        float pTowardEye;
        bool isSpecular;

    };
}

#endif // GOBLIN_PATH_VERTEX_H
