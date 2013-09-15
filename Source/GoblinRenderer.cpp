#include "GoblinRenderer.h"
#include "GoblinScene.h"
#include "GoblinRay.h"
#include "GoblinColor.h"
#include "GoblinCamera.h"
#include "GoblinFilm.h"
#include "GoblinSampler.h"

namespace Goblin {

    Renderer::Renderer():mSamples(NULL), mSampler(NULL) {}

    Renderer::~Renderer() {
        if(mSamples) {
            delete [] mSamples;
            mSamples = NULL;
        }
        if(mSampler) {
            delete [] mSampler;
            mSampler = NULL;
        }
    }

    void Renderer::render(Scene* scene) {
        const CameraPtr camera = scene->getCamera();
        // TODO sampler should instantiated based on the scene description
        // instead of this kind of combo getter + hardcode 1, 1
        Film* film = camera->getFilm();
        int xRes = film->getXResolution();
        int yRes = film->getYResolution();
        int xStart = film->getXStart();
        int yStart = film->getYStart();
        int xEnd = film->getXEnd();
        int yEnd = film->getYEnd();
        if(mSampler != NULL) {
            delete mSampler;
        }
        if(mSamples != NULL) {
            delete [] mSamples;
        }
        mSampler = new Sampler(xStart, xEnd, yStart, yEnd, 1, 1);
        mSamples = new CameraSample[mSampler->maxSamplesPerRequest()];
        int sampleNum = 0;
        while((sampleNum = mSampler->requestSamples(mSamples)) > 0) {
            for(int i = 0; i < sampleNum; ++i) {
                Ray ray;
                float w = camera->generateRay(mSamples[i], &ray);
                Color L = w * Li(scene, ray);
                film->addSample(mSamples[i], L);
            }
        }
    }

    Color Renderer::Li(Scene* scene, const Ray& ray) {
        Color Li = Color::Black;

        return Li;
    }
}
