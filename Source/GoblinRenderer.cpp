#include "GoblinRenderer.h"
#include "GoblinScene.h"
#include "GoblinRay.h"
#include "GoblinColor.h"
#include "GoblinCamera.h"
#include "GoblinFilm.h"
#include "GoblinSampler.h"

namespace Goblin {

    Renderer::Renderer(int maxRayDepth):
        mSamples(NULL), mSampler(NULL), mMaxRayDepth(maxRayDepth) {}

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
        //for now leave this here as debug purpose......
        //Ray ray(Vector3(0, 0, 0.0), Vector3(0, 0, 1), 0); 
        //Color L = Li(scene, ray);
        film->writeImage();
    }

    Color Renderer::Li(ScenePtr scene, const Ray& ray) const {
        Color Li = Color::Black;
        float epsilon;
        Intersection intersection;
        if(scene->intersect(ray, &epsilon, &intersection)) {
            const MaterialPtr& material = 
                intersection.primitive->getMaterial();
            const std::vector<LightPtr>& lights = scene->getLights();
            // direct light contribution
            for(size_t i = 0; i < lights.size(); ++i) {
                Vector3 wo = -ray.d;
                Vector3 wi;
                Ray shadowRay;
                const Fragment& fragment = intersection.fragment;
                Color L = lights[i]->Li(fragment.position, epsilon, &wi, 
                    &shadowRay);
                Color f = material->bsdf(fragment, wo, wi);
                if(f != Color::Black && !scene->intersect(shadowRay)) {
                    Li += f * L * absdot(fragment.normal, wi);
                }
            }
            // reflection and refraction
            if(ray.depth < mMaxRayDepth) {
                Li += specularReflect(scene, ray, epsilon, intersection);
                Li += specularRefract(scene, ray, epsilon, intersection);
            }
        }
        return Li;
    }

    Color Renderer::specularReflect(const ScenePtr& scene, const Ray& ray, 
        float epsilon, const Intersection& intersection) const {
        Color L(Color::Black);
        const Vector3& n = intersection.fragment.normal;
        const Vector3& p = intersection.fragment.position;
        const MaterialPtr& material = 
            intersection.primitive->getMaterial();
        Vector3 wo = -ray.d;
        Vector3 wi;
        Color f = material->sampleBSDF(intersection.fragment, 
            wo, &wi, BSDFType(BSDFSpecular | BSDFReflection));
        if(f != Color::Black && absdot(wi, n) != 0.0f) {
            Ray reflectiveRay(p, wi, epsilon);
            reflectiveRay.depth = ray.depth + 1;
            Color Lr = Li(scene, reflectiveRay);
            L += f * Lr * absdot(wi, n);
        }
        return L;
    }

    Color Renderer::specularRefract(const ScenePtr& scene, const Ray& ray, 
        float epsilon, const Intersection& intersection) const {
        Color L(Color::Black);
        const Vector3& n = intersection.fragment.normal;
        const Vector3& p = intersection.fragment.position;
        const MaterialPtr& material = 
            intersection.primitive->getMaterial();
        Vector3 wo = -ray.d;
        Vector3 wi;
        Color f = material->sampleBSDF(intersection.fragment, 
            wo, &wi, BSDFType(BSDFSpecular | BSDFTransmission));
        if(f != Color::Black && absdot(wi, n) != 0.0f) {
            Ray refractiveRay(p, wi, epsilon);
            refractiveRay.depth = ray.depth + 1;
            Color Lr = Li(scene, refractiveRay);
            L += f * Lr * absdot(wi, n);
        }
        return L;
    }

}
