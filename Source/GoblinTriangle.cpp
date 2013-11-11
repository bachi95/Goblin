#include "GoblinTriangle.h"
#include "GoblinObjMesh.h"
#include "GoblinRay.h"
#include "GoblinBBox.h"

namespace Goblin {
    Triangle::Triangle(ObjMesh* parentMesh, size_t index):
        mParentMesh(parentMesh), mIndex(index) {
        TriangleIndex* ti = (TriangleIndex*)parentMesh->getFacePtr(index);
        unsigned int i0 = ti->v[0];
        unsigned int i1 = ti->v[1];
        unsigned int i2 = ti->v[2];

        Vector3& v0 = ((Vertex*)parentMesh->getVertexPtr(i0))->position;
        Vector3& v1 = ((Vertex*)parentMesh->getVertexPtr(i1))->position;
        Vector3& v2 = ((Vertex*)parentMesh->getVertexPtr(i2))->position;
    }


    /*
     * solve ray.o + ray.d * t = some point in triangle
     * p(b1, b2) = (1 - b1 - b2) * v0 + b1 * v1 + b2 * v2
     * (barycentric coordinate b1 + b2 <=1 b1 >=0 b2 >= 0)
     * let e1 = v1 -v0, e2 = v2 - v0, s = ray.o - v0
     *                  |t |
     * |-ray.d e1 e2| * |b1| = s
     *                  |b2|
     * try to solve the above linear equation =>
     * according to Cramer's rule:
     * |t |                                | det|     s e1 e2| |
     * |b1| = (1.0f / det|-ray.d e1 e2|) * | det|-ray.d  s e2| |
     * |b2|                                | det|-ray.d e1  s| |
     * 
     * make use of the property:
     * det|a b c| = dot(a, cross(b, c)) = 
     *              dot(b, cross(c, a)) = 
     *              dot(c, cross(a, b))
     * let s1 = cross(ray.d, e2), s2 = cross(s, e1)
     * |t |                          | dot(s2,    e2) |
     * |b1| = (1.0f / dot(s1, e1)) * | dot(s1,     s) |
     * |b2|                          | dot(s2, ray.d) |
     * 
     * There is a more intuitive version of the above equation
     * basic idea first get intersection of ray and triangle plane, get t
     * then calculate b1, b2 with 2 2X2 linear equation, but the above 
     * one is cheaper to implement
     */
    bool Triangle::intersect(const Ray& ray) {

        TriangleIndex* ti = (TriangleIndex*)mParentMesh->getFacePtr(mIndex);
        unsigned int i0 = ti->v[0];
        unsigned int i1 = ti->v[1];
        unsigned int i2 = ti->v[2];
        Vector3& v0 = ((Vertex*)mParentMesh->getVertexPtr(i0))->position;
        Vector3& v1 = ((Vertex*)mParentMesh->getVertexPtr(i1))->position;
        Vector3& v2 = ((Vertex*)mParentMesh->getVertexPtr(i2))->position;

        Vector3 e1 = v1 - v0;
        Vector3 e2 = v2 - v0;
        Vector3 s1 = cross(ray.d, e2);
        float divisor = dot(s1, e1);
        if(divisor == 0.0f) {
            return false;
        }
        float invDivisor = 1.0f / divisor;

        Vector3 s = ray.o - v0;
        float b1 = dot(s, s1) * invDivisor;
        if(b1 < 0.0f || b1 > 1.0f) {
            return false;
        }

        Vector3 s2 = cross(s, e1);
        float b2 = dot(ray.d, s2) * invDivisor;
        if(b2 < 0.0f || b1 + b2 > 1.0f) {
            return false;
        }

        float t = dot(e2, s2) * invDivisor;
        if(t < ray.mint || t > ray.maxt) {
            return false;
        }

        return true;
    }

    bool Triangle::intersect(const Ray& ray, float* epsilon, 
        Intersection* intersection) {
        TriangleIndex* ti = (TriangleIndex*)mParentMesh->getFacePtr(mIndex);
        unsigned int i0 = ti->v[0];
        unsigned int i1 = ti->v[1];
        unsigned int i2 = ti->v[2];
        Vector3& v0 = ((Vertex*)mParentMesh->getVertexPtr(i0))->position;
        Vector3& v1 = ((Vertex*)mParentMesh->getVertexPtr(i1))->position;
        Vector3& v2 = ((Vertex*)mParentMesh->getVertexPtr(i2))->position;

        Vector3 e1 = v1 - v0;
        Vector3 e2 = v2 - v0;
        Vector3 s1 = cross(ray.d, e2);
        float divisor = dot(s1, e1);
        if(divisor == 0.0f) {
            return false;
        }
        float invDivisor = 1.0f / divisor;

        Vector3 s = ray.o - v0;
        float b1 = dot(s, s1) * invDivisor;
        if(b1 < 0.0f || b1 > 1.0f) {
            return false;
        }

        Vector3 s2 = cross(s, e1);
        float b2 = dot(ray.d, s2) * invDivisor;
        if(b2 < 0.0f || b1 + b2 > 1.0f) {
            return false;
        }

        float t = dot(e2, s2) * invDivisor;
        if(t < ray.mint || t > ray.maxt) {
            return false;
        }

        ray.maxt = t;
        *epsilon = 1e-3f * t;
        intersection->position = ray(t);
        Vector3& n0 = ((Vertex*)mParentMesh->getVertexPtr(i0))->normal;
        Vector3& n1 = ((Vertex*)mParentMesh->getVertexPtr(i1))->normal;
        Vector3& n2 = ((Vertex*)mParentMesh->getVertexPtr(i2))->normal;
        intersection->normal = (1.0f - b1 - b2) * n0 + b1 * n1 + b2 * n2;
        intersection->normal.normalize();
        return true;
    }

    BBox Triangle::getObjectBound() {
        TriangleIndex* ti = (TriangleIndex*)mParentMesh->getFacePtr(mIndex);
        unsigned int i0 = ti->v[0];
        unsigned int i1 = ti->v[1];
        unsigned int i2 = ti->v[2];
        Vector3& v0 = ((Vertex*)mParentMesh->getVertexPtr(i0))->position;
        Vector3& v1 = ((Vertex*)mParentMesh->getVertexPtr(i1))->position;
        Vector3& v2 = ((Vertex*)mParentMesh->getVertexPtr(i2))->position;

        BBox rv;
        rv.expand(v0);
        rv.expand(v1);
        rv.expand(v2);
        return rv;
    }
}
