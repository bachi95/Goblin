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

    void Renderer::render(ScenePtr scene) {
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
        film->writeImage();
    }

    Color Renderer::Li(ScenePtr scene, const Ray& ray) {
        Color Li = Color::Black;
        float epsilon;
        Intersection intersection;
        if(scene->intersect(ray, &epsilon, &intersection)) {
            const MaterialPtr& material = 
                intersection.primitive->getMaterial();
            const std::vector<LightPtr>& lights = scene->getLights();
            for(size_t i = 0; i < lights.size(); ++i) {
                Vector3 wi;
                Ray shadowRay;
                Color L = lights[i]->Li(intersection.position, epsilon, &wi, 
                    &shadowRay);
                Color f = material->bsdf(Vertex(), wi, -ray.d);
                if(f != Color::Black && !scene->intersect(shadowRay)) {
                    Li += f * L * clamp(dot(intersection.normal, wi), 0, 1);
                }
            }
        }
        return Li;
    }
}
