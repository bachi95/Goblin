#include "GoblinSphere.h"
#include "GoblinUtils.h"
#include "GoblinRay.h"
#include "GoblinBBox.h"
#include "GoblinSampler.h"

using namespace Goblin;

Sphere::Sphere(float r, size_t numSlices, size_t numStacks):
    mRadius(r), 
    mNumSlices(numSlices), 
    mNumStacks(numStacks) { }

void Sphere::init() {
    buildStacks();
}

bool Sphere::intersect(const Ray& ray) const {
    float A = squaredLength(ray.d);
    float B = 2.0f * dot(ray.d, ray.o);
    float C = squaredLength(ray.o) - mRadius * mRadius;
    float tNear, tFar;
    if(!quadratic(A, B, C, &tNear, &tFar)) {
        return false;
    }
    if(tNear > ray.maxt || tFar < ray.mint) {
        return false;
    }
    float tHit = tNear;
    if(tHit < ray.mint) {
        tHit = tFar;
        if(tHit > ray.maxt) {
            return false;
        }
    }
    return true;
}

bool Sphere::intersect(const Ray& ray, float* epsilon, 
        Fragment* fragment) const {
    float A = squaredLength(ray.d);
    float B = 2.0f * dot(ray.d, ray.o);
    float C = squaredLength(ray.o) - mRadius * mRadius;
    float tNear, tFar;
    if(!quadratic(A, B, C, &tNear, &tFar)) {
        return false;
    }
    if(tNear > ray.maxt || tFar < ray.mint) {
        return false;
    }
    float tHit = tNear;
    if(tHit < ray.mint) {
        tHit = tFar;
        if(tHit > ray.maxt) {
            return false;
        }
    }
    ray.maxt = tHit;
    *epsilon = 1e-3f * tHit;
    Vector3 pHit = ray(tHit);
    fragment->position = pHit;
    fragment->normal = normalize(fragment->position);

    /* math time~
    spherical cooridnate
    x = r * sinTheta * cosPhi
    y = r * cosTheta
    z = r * sinTheta * sinPhi
    0 < = Phi <= 2PI, 0 <= Theta <=PI
    u = Phi / 2PI
    v = Theta / PI

    dPx/du = d(r * sinTheta * cosPhi)/du = r * sinTheta * d(cos(2PI * u))/du =
    r * sinTheta *-sin(2PI * u) * 2PI = -2PI * z
    dPx/du = d(r * cosTheta)/du = 0
    dPx/dz = d(r * sinTheta * sinPhi)/du = r * sinTheta * d(sin(2PI * u))/du =
    r * sinTheta *cos(2PI * u) * 2PI = 2PI * x
    dpdu = [-2PI * z, 0, 2PI * x]

    dPx/dv = d(r * sinTheta * cosPhi)/dv = r * cosPhi * d(sin(PI * v))/dv =
    r * cosPhi * cos(PI * v) * PI = y * PI * cosPhi
    dPy/dv = d(r * cosTheta)/dv = r * d(cos(PI * v))/dv = r * PI * -sinTheta =
    -PI * r * sinTheta
    dPz/dv = d(r * sinTheta * sinPhi)/dv = r * sinPhi * d(sin(PI * v))/dv =
    r * sinPhi * cos(PI * v) * PI = y * PI * sinPhi
    dpdv = PI * (y * cosPhi, -r * sinTheta, y * sinPhi)

    cosPhi = x / sqrt(x * x + z * z)
    sinPhi = z / sqrt(x * x + z * z)
    */
    float phi = atan2(pHit.z, pHit.x); 
    if(phi < 0.0f) {
        phi += 2.0f * PI;
    }
    float u = phi * INV_TWOPI;
    float theta = acos(pHit.y / mRadius);
    float v = theta * INV_PI;
    float invR = 1.0f / sqrt(pHit.x * pHit.x + pHit.z * pHit.z);
    float cosPhi = pHit.x * invR;
    float sinPhi = pHit.z * invR;
    fragment->uv = Vector2(u, v);
    fragment->dpdu = Vector3(-TWO_PI * pHit.z, 0.0f, TWO_PI * pHit.x);
    fragment->dpdv = PI * Vector3(pHit.y * cosPhi, -mRadius * sin(theta),
        pHit.y * sinPhi);
    return true;
}

Vector3 Sphere::sample(float u1, float u2, Vector3* normal) const {
    *normal = mRadius * uniformSampleSphere(u1, u2);
    return mRadius * (*normal);
}

BBox Sphere::getObjectBound() {
    return BBox(Vector3(mRadius, mRadius, mRadius), 
        Vector3(-mRadius, -mRadius, -mRadius));
}

void Sphere::buildStacks() {
    float thetaStep = PI / mNumStacks;
    float phiStep = 2.0f * PI / mNumSlices;
    // two poles of the sphere are not counted as ring
    size_t numRings = mNumStacks - 1;
    for(size_t i = 1; i <= numRings; ++i) {
        float theta = i * thetaStep;
        for(size_t j = 0; j <= mNumSlices; ++j) {
            float phi = j * phiStep;
            // from top to bottom
            Vertex v;
            float sinTheta = sin(theta);
            float cosTheta = cos(theta);
            float sinPhi = sin(phi);
            float cosPhi = cos(phi);
            v.position.x = mRadius * sinTheta * cosPhi;
            v.position.y = mRadius * cosTheta;
            v.position.z = mRadius * sinTheta * sinPhi;

            v.tangent.x = -mRadius * cosTheta * sinPhi;
            v.tangent.y = 0.0f;
            v.tangent.z = mRadius * sinTheta * cosPhi; 
            
            v.normal = normalize(v.position);

            v.texC.x = phi * INV_TWOPI;
            v.texC.y = theta * INV_PI;

            mVertices.push_back(v);
        }
    }

    mVertices.push_back(Vertex(0.0f, -mRadius, 0.0f, 
        1.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 0.0f, 1.0f));
    mVertices.push_back(Vertex(0.0f, mRadius, 0.0f,
        1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f));
    size_t northPoleIndex = mVertices.size() - 1;
    size_t southPoleIndex = mVertices.size() - 2;

    size_t numRingVertices = mNumSlices + 1;

    for(size_t i = 0; i < mNumStacks - 2; ++i) {
        for(size_t j = 0; j < mNumSlices; ++j) {
            TriangleIndex triangle;
            triangle.v[0] = i * numRingVertices + j;
            triangle.v[1] = i * numRingVertices + j + 1;
            triangle.v[2] = (i + 1) * numRingVertices + j;
            mTriangles.push_back(triangle);

            triangle.v[0] = (i + 1) * numRingVertices + j;
            triangle.v[1] = i * numRingVertices + j + 1;
            triangle.v[2] = (i + 1) * numRingVertices + j + 1;
            mTriangles.push_back(triangle);
        }
    }

    //top ring is in the rear block of vertices
    for(size_t i = 0; i < mNumSlices; ++i) {
        TriangleIndex triangle;
        triangle.v[0] = northPoleIndex;
        triangle.v[1] = i + 1;
        triangle.v[2] = i;
        mTriangles.push_back(triangle);
    }
    // how the vertices layout:
    // | top ring . middle rings . bottom ring. south pole, north pole |
    size_t baseIndex = (numRings - 1) * numRingVertices;
    for(size_t i = 0; i < mNumSlices; ++i) {
        TriangleIndex triangle;
        triangle.v[0] = southPoleIndex;
        triangle.v[1] = baseIndex + i;
        triangle.v[2] = baseIndex + i + 1;
        mTriangles.push_back(triangle);
    }
}
