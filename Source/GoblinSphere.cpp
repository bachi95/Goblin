#include "GoblinSphere.h"
#include "GoblinUtils.h"
#include "GoblinRay.h"

using namespace Goblin;

Sphere::Sphere(float r, size_t numSlices, size_t numStacks):
    mRadius(r), 
    mNumSlices(numSlices), 
    mNumStacks(numStacks) { }

void Sphere::init() {
    buildStacks();
}

bool Sphere::intersect(const Ray& ray) {
    float A = squaredLength(ray.d);
    float B = 2.0f * dot(ray.d, ray.o);
    float C = squaredLength(ray.o) - mRadius * mRadius;
    float tNear, tFar;
    if(!quadratic(A, B, C, &tNear, &tFar)) {
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

void Sphere::buildStacks() {
    float phiStep = PI / mNumStacks;
    float thetaStep = 2.0f * PI / mNumSlices;
    // two poles of the sphere are not counted as ring
    size_t numRings = mNumStacks - 1;
    for(size_t i = 1; i <= numRings; ++i) {
        float phi = i * phiStep;
        for(size_t j = 0; j <= mNumSlices; ++j) {
            float theta = j * thetaStep;
            // from top to bottom
            Vertex v;
            v.position.x = mRadius * sinf(phi) * cosf(theta);
            v.position.y = mRadius * cosf(phi);
            v.position.z = mRadius * sinf(phi) * sinf(theta);

            v.tangent.x = -mRadius * sinf(phi) * sinf(theta);
            v.tangent.y = 0.0f;
            v.tangent.z = mRadius * sinf(phi) * cosf(theta); 
            
            v.normal = normalize(v.position);

            v.texC.x = theta / (2.0f * PI);
            v.texC.y = phi / PI;

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
