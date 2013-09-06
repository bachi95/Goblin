#include "GoblinRenderer.h"
#include "GoblinScene.h"
#include "GoblinRay.h"
#include "GoblinColor.h"
#include "GoblinCamera.h"
#include "GoblinFilm.h"
#include "GoblinSampler.h"

namespace Goblin {
    void Renderer::render(Scene* scene) {
        const CameraPtr camera = scene->getCamera();
        //TODO replace this naive X*Y ray shooting with sampler
        Film* film = camera->getFilm();
        int xRes = film->getXResolution();
        int yRes = film->getYResolution();
        for(int i = 0; i < yRes; ++i) {
            for(int j = 0; j < xRes; ++j) {
                CameraSample sample((float)i, (float)j);
                Ray ray;
                float w = camera->generateRay(sample, &ray);
                Color L = w * Li(scene, ray);
                film->addSample(sample, L);
            }
        }
    }

    Color Renderer::Li(Scene* scene, const Ray& ray) {
        Color Li = Color::Black;

        return Li;
    }
}