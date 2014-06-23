#include "GoblinImageIO.h"
#include "GoblinLight.h"
#include "GoblinRay.h"
#include "GoblinSampler.h"
#include "GoblinScene.h"

namespace Goblin {

    LightSampleIndex::LightSampleIndex(SampleQuota* sampleQuota, 
        int requestNum) {
        SampleIndex oneDIndex = sampleQuota->requestOneDQuota(requestNum);
        SampleIndex twoDIndex = sampleQuota->requestTwoDQuota(requestNum);
        // theoretically this two should be the same...
        // this is just a paranoid double check
        samplesNum = min(oneDIndex.sampleNum, twoDIndex.sampleNum);
        componentIndex = oneDIndex.offset;
        geometryIndex = twoDIndex.offset;
    }

    LightSample::LightSample(const RNG& rng) {
        uComponent = rng.randomFloat();
        uGeometry[0] = rng.randomFloat();
        uGeometry[1] = rng.randomFloat();
    }

    LightSample::LightSample(const Sample& sample, 
        const LightSampleIndex& index, uint32_t n) {
        uComponent = sample.u1D[index.componentIndex][n];
        uGeometry[0] = sample.u2D[index.geometryIndex][2 * n];
        uGeometry[1] = sample.u2D[index.geometryIndex][2 * n + 1];
    }

    PointLight::PointLight(const Color& I, const Vector3& P):
    intensity(I), position(P) {
        mParams.setInt("type", Point);
        mParams.setColor("intensity", intensity);
        mParams.setVector3("position", position);
    }

    Color PointLight::sampleL(const Vector3& p, float epsilon, 
        const LightSample& lightSample,
        Vector3* wi, float* pdf, Ray* shadowRay) const {
        Vector3 dir = position - p;
        *wi = normalize(dir);
        *pdf = 1.0f;
        shadowRay->o = p;
        shadowRay->d = *wi;
        shadowRay->mint = epsilon;
        float squaredDistance = squaredLength(dir);
        shadowRay->maxt = sqrt(squaredDistance) - epsilon;
        return intensity / squaredDistance;
    }

    Color PointLight::sampleL(const ScenePtr& scene, const LightSample& ls,
        float u1, float u2, Ray* ray, float* pdf) const {
        Vector3 dir = uniformSampleSphere(ls.uGeometry[0], ls.uGeometry[1]);
        *ray = Ray(position, dir, 0.0f);
        if(pdf) {
            *pdf = uniformSpherePdf();
        }
        return intensity;
    }

    Color PointLight::power(const ScenePtr& scene) const {
        return 4.0f * PI * intensity;
    }

    DirectionalLight::DirectionalLight(const Color& R, const Vector3& D):
    radiance(R), direction(D) {
        mParams.setInt("type", Directional);
        mParams.setColor("radiance", radiance);
        mParams.setVector3("direction", direction);
    }

    Color DirectionalLight::sampleL(const Vector3& p, float epsilon, 
        const LightSample& lightSample,
        Vector3* wi, float* pdf, Ray* shadowRay) const {
        *wi = -direction;
        *pdf = 1.0f;
        shadowRay->o = p;
        shadowRay->d = *wi;
        shadowRay->mint = epsilon;
        return radiance;
    }

    /*
     * approximation of sample directional light by sampling among
     * the world bounding sphere, first sample a point from disk
     * with world radius that perpendicular to light direction, 
     * then offset it back world radius distance as ray origin
     * ray dir is simply light dir
     */
    Color DirectionalLight::sampleL(const ScenePtr& scene, 
        const LightSample& ls,
        float u1, float u2, Ray* ray, float* pdf) const {
        Vector3 worldCenter;
        float worldRadius;
        scene->getBoundingSphere(&worldCenter, &worldRadius);
        Vector3 xAxis, yAxis;
        coordinateAxises(direction, &xAxis, &yAxis);
        Vector2 diskXY = uniformSampleDisk(ls.uGeometry[0], ls.uGeometry[1]);
        Vector3 worldDiskSample = worldCenter + 
            worldRadius * (diskXY.x * xAxis + diskXY.y * yAxis);
        Vector3 origin = worldDiskSample - direction * worldRadius;
        *ray = Ray(origin, direction, 0.0f);
        if(pdf) {
            *pdf = INV_PI / (worldRadius * worldRadius);
        }
        return radiance;
    }

    Color DirectionalLight::power(const ScenePtr& scene) const {
        Vector3 center;
        float radius;
        scene->getBoundingSphere(&center, &radius);
        // well...... we can't make it infinitely big, so use bounding
        // sphere for a rough approximation 
        return radius * radius * PI * radiance;
    }

    GeometrySet::GeometrySet(const GeometryPtr& geometry):
        mSumArea(0.0f), mAreaDistribution(NULL) {
        if(geometry->intersectable()) {
            mGeometries.push_back(geometry);
        } else {
            geometry->refine(mGeometries);
        }
        mSumArea = 0.0f;
        mGeometriesArea.resize(mGeometries.size());
        for(size_t i = 0; i < mGeometries.size(); ++i) {
            float area = mGeometries[i]->area();
            mGeometriesArea[i] = area;
            mSumArea += area;
        }
        mAreaDistribution = new CDF1D(mGeometriesArea);
    }

    GeometrySet::~GeometrySet() {
        if(mAreaDistribution != NULL ) {
            delete mAreaDistribution;
            mAreaDistribution = NULL;
        }
    }

    Vector3 GeometrySet::sample(const Vector3& p, 
        const LightSample& lightSample,
        Vector3* normal) const {
        // pick up a geometry to sample based on area distribution
        float uComp = lightSample.uComponent;
        int geoIndex = mAreaDistribution->sampleDiscrete(uComp);
        // sample out ps from picked up geometry surface
        float u1 = lightSample.uGeometry[0];
        float u2 = lightSample.uGeometry[1];
        Vector3 ps = mGeometries[geoIndex]->sample(p, u1, u2, normal);
        return ps;
    }

    Vector3 GeometrySet::sample(const LightSample& lightSample,
        Vector3* normal) const {
        float uComp = lightSample.uComponent;
        int geoIndex = mAreaDistribution->sampleDiscrete(uComp);
        float u1 = lightSample.uGeometry[0];
        float u2 = lightSample.uGeometry[1];
        Vector3 ps = mGeometries[geoIndex]->sample(u1, u2, normal);
        return ps;
    }

    float GeometrySet::pdf(const Vector3& p, const Vector3& wi) const {
        float pdf = 0.0f;
        for(size_t i = 0; i < mGeometries.size(); ++i) {
            pdf += mGeometriesArea[i] * mGeometries[i]->pdf(p, wi);    
        }
        pdf /= mSumArea;
        return pdf;
    }


    AreaLight::AreaLight(const Color& Le, const GeometryPtr& geometry,
        const Transform& toWorld, uint32_t samplesNum): mLe(Le), 
        mToWorld(toWorld), mSamplesNum(samplesNum) {
        mGeometrySet = new GeometrySet(geometry);
    }

    AreaLight::~AreaLight() {
        if(mGeometrySet != NULL) {
            delete mGeometrySet;
            mGeometrySet = NULL;
        }
    }

    Color AreaLight::L(const Vector3& ps, const Vector3& ns, 
        const Vector3& w) const {
        return dot(ns, w) > 0.0f ? mLe : Color::Black;
    }

    Color AreaLight::sampleL(const Vector3& p, float epsilon,
        const LightSample& lightSample,
        Vector3* wi, float* pdf, Ray* shadowRay) const {
        // transform world space p to local space since all GeometrySet methods
        // are in local space
        Vector3 pLocal = mToWorld.invertPoint(p);
        Vector3 nsLocal;
        Vector3 psLocal = mGeometrySet->sample(pLocal, lightSample, &nsLocal);
        Vector3 wiLocal = normalize(psLocal - pLocal);
        *pdf = mGeometrySet->pdf(pLocal, wiLocal);
        // transform
        Vector3 ps = mToWorld.onPoint(psLocal); 
        Vector3 ns = normalize(mToWorld.onNormal(nsLocal));
        *wi = normalize(ps - p);

        shadowRay->o = p;
        shadowRay->d = *wi;
        shadowRay->mint = epsilon;
        shadowRay->maxt = length(ps - p) - epsilon;

        return L(ps, ns, -*wi);
    }

    Color AreaLight::sampleL(const ScenePtr& scene, const LightSample& ls,
        float u1, float u2, Ray* ray, float* pdf) const {
        Vector3 n;
        Vector3 origin = mGeometrySet->sample(ls, &n);
        Vector3 dir = uniformSampleSphere(u1, u2);
        // the case sampled dir is in the opposite hemisphere of area light
        // surface normal
        if(dot(n, dir) < 0.0f) {
            dir *= -1.0f;
        }
        *ray = Ray(origin, dir, 1e-3f);
        // each point on the area light has uniform hemisphere distribution
        // output radiance
        if(pdf) {
            Vector3 worldScale = mToWorld.getScale();
            float worldArea = mGeometrySet->area() *
                worldScale.x * worldScale.y * worldScale.z;
            *pdf = (1.0f / worldArea) * INV_TWOPI;
        }
        return mLe;
    }

    Color AreaLight::power(const ScenePtr& scene) const {
        // if any random angle output mLe on the area light surface,
        // we can think of the input radience with perpendular angle
        // per unit area mLe * PI (similar to how we get lambert bsdf)
        return mLe * PI * mGeometrySet->area();
    }

    float AreaLight::pdf(const Vector3& p, const Vector3& wi) const {
        Vector3 pLocal = mToWorld.invertPoint(p);
        Vector3 wiLocal = mToWorld.invertVector(wi);
        return mGeometrySet->pdf(pLocal, wiLocal);
    }


    ImageBasedLight::ImageBasedLight(const string& radianceMap, 
        const Color& filter, const Quaternion& orientation): 
        mRadiance(NULL), mDistribution(NULL) {
        // Make default orientation facing the center of
        // environment map since spherical coordinate is z-up
        mToWorld.rotateX(-0.5f * PI);
        mToWorld.rotateY(-0.5f * PI);
        mToWorld.setOrientation(orientation * mToWorld.getOrientation());
        int width, height;
        Color* buffer = loadImage(radianceMap, &width, &height);
        if(buffer == NULL) {
            std::cerr << "errror loading image " << radianceMap << std::endl;
            width = 1;
            height = 1;
            buffer = new Color[1];
            buffer[0] = Color::Magenta;
        }
        for(int i = 0; i < width * height; ++i) {
            buffer[i] *= filter;
            mAverageRadiance += buffer[i];
        }
        mAverageRadiance /= (float)(width * height);
        mRadiance = new ImageBuffer<Color>(buffer, width, height);
        float* dist = new float[width * height];
        for(int i = 0; i < height; ++i) {
            float sinTheta = sin(((float)i + 0.5f) / (float)height * PI);
            for(int j = 0; j < width; ++j) {
                int index = i * width + j;
                dist[index] = buffer[index].luminance() * sinTheta;
            }
        }
        mDistribution = new CDF2D(dist, width, height);
        delete [] dist;
    }

    ImageBasedLight::~ImageBasedLight() {
        if(mRadiance != NULL) {
            delete mDistribution;
            mDistribution = NULL;
        }
        if(mDistribution != NULL) {
            delete mDistribution;
            mDistribution = NULL;
        }
    }

    Color ImageBasedLight::Le(const Ray& ray) const {
        const Vector3& w = mToWorld.invertVector(ray.d);
        float theta = acos(w.z);
        float phi = atan2(w.y, w.x);
        if(phi < 0.0f) {
            phi += TWO_PI;
        }
        float s = phi * INV_TWOPI;
        float t = theta * INV_PI;
        return mRadiance->lookup(s, t);
    }

    Color ImageBasedLight::sampleL(const Vector3& p, float epsilon,
        const LightSample& lightSample,
        Vector3* wi, float* pdf, Ray* shadowRay) const {
        float pdfST;
        Vector2 st = mDistribution->sampleContinuous(
            lightSample.uGeometry[0], lightSample.uGeometry[1], &pdfST);
        float theta = st[1] * PI;
        float phi = st[0] * TWO_PI;
        float cosTheta = cos(theta);
        float sinTheta = sin(theta);
        float cosPhi = cos(phi);
        float sinPhi = sin(phi);
        Vector3 wLocal(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);
        *wi = mToWorld.onVector(wLocal);

        if(sinTheta == 0.0f) {
            *pdf = 0.0f;
        }
        *pdf = pdfST / (TWO_PI * PI * sinTheta);

        shadowRay->o = p;
        shadowRay->d = *wi;
        shadowRay->mint = epsilon;

        return mRadiance->lookup(st[0], st[1]);
    }

    Color ImageBasedLight::sampleL(const ScenePtr& scene, 
        const LightSample& ls, float u1, float u2, 
        Ray* ray, float* pdf) const {

        float pdfST;
        Vector2 st = mDistribution->sampleContinuous(
            ls.uGeometry[0], ls.uGeometry[1], &pdfST);
        float theta = st[1] * PI;
        float phi = st[0] * TWO_PI;
        float cosTheta = cos(theta);
        float sinTheta = sin(theta);
        float cosPhi = cos(phi);
        float sinPhi = sin(phi);
        Vector3 w = mToWorld.onVector(
            Vector3(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta));

        Vector3 worldCenter;
        float worldRadius;
        scene->getBoundingSphere(&worldCenter, &worldRadius);
        Vector3 xAxis, yAxis;
        coordinateAxises(-w, &xAxis, &yAxis);
        Vector2 diskXY = uniformSampleDisk(u1, u2);
        Vector3 worldDiskSample = worldCenter + 
            worldRadius * (diskXY.x * xAxis + diskXY.y * yAxis);
        Vector3 origin = worldRadius * w;
        Vector3 direction = normalize(worldDiskSample - origin);
        *ray = Ray(origin, direction, 0.0f);
        if(pdf) {
            *pdf = pdfST / (TWO_PI * PI * sinTheta);
        }
        return mRadiance->lookup(st[0], st[1]);
    }

    Color ImageBasedLight::power(const ScenePtr& scene) const {
        Vector3 center;
        float radius;
        scene->getBoundingSphere(&center, &radius);
        // raough power estimation, assume radiance in world sphere
        // diffuse distribution
        return mAverageRadiance * PI * (4.0f * PI * radius * radius);
    }

    float ImageBasedLight::pdf(const Vector3& p, const Vector3& wi) const {
        Vector3 wiLocal = mToWorld.invertVector(wi);
        float theta = acos(wiLocal.z);
        float sinTheta = sin(theta);
        if(sinTheta == 0.0f) {
            return 0.0f;
        }

        float phi = atan2(wi.y, wi.x);
        if(phi < 0.0f) {
            phi += TWO_PI;
        }
        float pdf = mDistribution->pdf(phi * INV_TWOPI, theta * INV_PI) / 
            (TWO_PI * PI * sinTheta);
        return pdf;
    }
}
