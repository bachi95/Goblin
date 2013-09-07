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
        for(int y = 0; y < yRes; ++y) {
            for(int x = 0; x < xRes; ++x) {
                CameraSample sample((float)x, (float)y);
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
