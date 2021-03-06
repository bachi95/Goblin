#include "GoblinTriangle.h"
#include "GoblinPolygonMesh.h"
#include "GoblinRay.h"
#include "GoblinBBox.h"
#include "GoblinSampler.h"

namespace Goblin {
Triangle::Triangle(const PolygonMesh* parentMesh, size_t index):
    mParentMesh(parentMesh), mIndex(index) {}


// solve ray.o + ray.d * t = some point in triangle
// p(b1, b2) = (1 - b1 - b2) * p0 + b1 * p1 + b2 * p2
// (barycentric coordinate b1 + b2 <=1 b1 >=0 b2 >= 0)
// let e1 = p1 -p0, e2 = p2 - p0, s = ray.o - p0
//                  |t |
// |-ray.d e1 e2| * |b1| = s
//                  |b2|
// try to solve the above linear equation =>
// according to Cramer's rule:
// |t |                                | det|     s e1 e2| |
// |b1| = (1.0f / det|-ray.d e1 e2|) * | det|-ray.d  s e2| |
// |b2|                                | det|-ray.d e1  s| |
// 
// make use of the property:
// det|a b c| = dot(a, cross(b, c)) = 
//             dot(b, cross(c, a)) = 
//             dot(c, cross(a, b))
// let s1 = cross(ray.d, e2), s2 = cross(s, e1)
// |t |                          | dot(s2,    e2) |
// |b1| = (1.0f / dot(s1, e1)) * | dot(s1,     s) |
// |b2|                          | dot(s2, ray.d) |
// 
// There is a more intuitive version of the above equation
// basic idea first get intersection of ray and triangle plane, get t
// then calculate b1, b2 with 2 2X2 linear equation, but the above 
// one is cheaper to implement
bool Triangle::intersect(const Ray& ray, float* epsilon, 
    Fragment* fragment) const {
    TriangleIndex* ti = (TriangleIndex*)mParentMesh->getFacePtr(mIndex);
    unsigned int i0 = ti->v[0];
    unsigned int i1 = ti->v[1];
    unsigned int i2 = ti->v[2];

    const Vertex& v0 = *mParentMesh->getVertexPtr(i0);
    const Vertex& v1 = *mParentMesh->getVertexPtr(i1);
    const Vertex& v2 = *mParentMesh->getVertexPtr(i2);

    const Vector3& p0 = v0.position;
    const Vector3& p1 = v1.position;
    const Vector3& p2 = v2.position;

    Vector3 e1 = p1 - p0;
    Vector3 e2 = p2 - p0;
    Vector3 s1 = cross(ray.d, e2);
    float divisor = dot(s1, e1);
    if (divisor == 0.0f) {
        return false;
    }
    float invDivisor = 1.0f / divisor;
    float fEpsilon = 1e-7f;
    Vector3 s = ray.o - p0;
    float b1 = dot(s, s1) * invDivisor;
    if (b1 + fEpsilon < 0.0f || b1 - fEpsilon > 1.0f) {
        return false;
    }

    Vector3 s2 = cross(s, e1);
    float b2 = dot(ray.d, s2) * invDivisor;
    if (b2 + fEpsilon < 0.0f || b1 + b2 - fEpsilon > 1.0f) {
        return false;
    }

    float t = dot(e2, s2) * invDivisor;
    if (t < ray.mint || t > ray.maxt) {
        return false;
    }

    float b0 = 1.0f - b1 - b2;
    ray.maxt = t;
    *epsilon = 1e-3f * t;
    // start collect intersection geometry info:
    // position, normal, uv, dpdu, dpdv....
    Vector3 position(ray(t));
    Vector3 normal;
    if (mParentMesh->hasNormal()) {
        const Vector3& n0 = v0.normal;
        const Vector3& n1 = v1.normal;
        const Vector3& n2 = v2.normal;
        normal = normalize(b0 * n0 + b1 * n1 + b2 * n2);
    } else {
        normal = normalize(cross(e1, e2));
    }

    Vector2 uvs[3];
    if (mParentMesh->hasTexCoord()) {
        uvs[0] = v0.texC;
        uvs[1] = v1.texC;
        uvs[2] = v2.texC;
    } else {
        uvs[0] = Vector2(0.0f, 0.0f);
        uvs[1] = Vector2(1.0f, 0.0f);
        uvs[2] = Vector2(0.0f, 1.0f);
    }
    Vector2 uv(b0 * uvs[0] + b1 * uvs[1] + b2 * uvs[2]);

    float du1 = uvs[1].x - uvs[0].x;
    float dv1 = uvs[1].y - uvs[0].y;
    float du2 = uvs[2].x - uvs[0].x;
    float dv2 = uvs[2].y - uvs[0].y;
    float determinant = du1 * dv2 - dv1 * du2;
    Vector3 dpdu, dpdv;
    if (determinant == 0.0f) {
        // form a random shading coordinate from normal then
        dpdu = normalize(e1 - 
            dot(fragment->getNormal(), e1) * fragment->getNormal());
        dpdv = cross(fragment->getNormal(), fragment->getDPDV());
    } else {
        float invDet = 1.0f / determinant;
        dpdu = invDet * (dv2 * e1 - dv1 * e2);
        dpdv = invDet * (-du2 * e1 + du1 * e2);
    }
    *fragment = Fragment(position, normal, uv, dpdu, dpdv);
    return true;
}

bool Triangle::occluded(const Ray& ray) const {

	TriangleIndex* ti = (TriangleIndex*)mParentMesh->getFacePtr(mIndex);
	unsigned int i0 = ti->v[0];
	unsigned int i1 = ti->v[1];
	unsigned int i2 = ti->v[2];
	const Vector3& p0 = mParentMesh->getVertexPtr(i0)->position;
	const Vector3& p1 = mParentMesh->getVertexPtr(i1)->position;
	const Vector3& p2 = mParentMesh->getVertexPtr(i2)->position;

	Vector3 e1 = p1 - p0;
	Vector3 e2 = p2 - p0;
	Vector3 s1 = cross(ray.d, e2);
	float divisor = dot(s1, e1);
	if (divisor == 0.0f) {
		return false;
	}
	float invDivisor = 1.0f / divisor;
	float fEpsilon = 1e-7f;
	Vector3 s = ray.o - p0;
	float b1 = dot(s, s1) * invDivisor;
	if (b1 + fEpsilon < 0.0f || b1 - fEpsilon > 1.0f) {
		return false;
	}

	Vector3 s2 = cross(s, e1);
	float b2 = dot(ray.d, s2) * invDivisor;
	if (b2 + fEpsilon < 0.0f || b1 + b2 - fEpsilon > 1.0f) {
		return false;
	}

	float t = dot(e2, s2) * invDivisor;
	if (t < ray.mint || t > ray.maxt) {
		return false;
	}
	return true;
}

Vector3 Triangle::sample(float u1, float u2, Vector3* normal) const {
    float b0, b1;
    uniformSampleTriangle(u1, u2, &b0, &b1);
    TriangleIndex* ti = (TriangleIndex*)mParentMesh->getFacePtr(mIndex);
    unsigned int i0 = ti->v[0];
    unsigned int i1 = ti->v[1];
    unsigned int i2 = ti->v[2];
    const Vector3& p0 = mParentMesh->getVertexPtr(i0)->position;
    const Vector3& p1 = mParentMesh->getVertexPtr(i1)->position;
    const Vector3& p2 = mParentMesh->getVertexPtr(i2)->position;
    *normal = normalize(cross(p1 - p0, p2 - p0));
    return b0 * p0 + b1 * p1 + (1.0f - b0 - b1) * p2;
}

inline float Triangle::area() const {
    TriangleIndex* ti = (TriangleIndex*)mParentMesh->getFacePtr(mIndex);
    unsigned int i0 = ti->v[0];
    unsigned int i1 = ti->v[1];
    unsigned int i2 = ti->v[2];
    const Vector3& p0 = mParentMesh->getVertexPtr(i0)->position;
    const Vector3& p1 = mParentMesh->getVertexPtr(i1)->position;
    const Vector3& p2 = mParentMesh->getVertexPtr(i2)->position;

    return 0.5f * length(cross(p1 - p0, p2 - p0));
}

BBox Triangle::getObjectBound() const {
    TriangleIndex* ti = (TriangleIndex*)mParentMesh->getFacePtr(mIndex);
    unsigned int i0 = ti->v[0];
    unsigned int i1 = ti->v[1];
    unsigned int i2 = ti->v[2];
    const Vector3& v0 = mParentMesh->getVertexPtr(i0)->position;
    const Vector3& v1 = mParentMesh->getVertexPtr(i1)->position;
    const Vector3& v2 = mParentMesh->getVertexPtr(i2)->position;

    BBox rv;
    rv.expand(v0);
    rv.expand(v1);
    rv.expand(v2);
    return rv;
}

}