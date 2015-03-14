#include "GoblinGeometry.h"
#include "GoblinTransform.h"
#include "GoblinRay.h"

namespace Goblin {
    std::map<size_t, Geometry*> Geometry::geometryCache;

    Fragment::Fragment(const Vector3& p, const Vector3& n, const Vector2& uv,
        const Vector3& dpdu, const Vector3& dpdv):
        mPosition(p), mNormal(n), mUV(uv), mDPDU(dpdu), mDPDV(dpdv),
        mDPDX(Vector3::Zero), mDPDY(Vector3::Zero),
        mDUDX(0.0f), mDVDX(0.0f), mDUDY(0.0f), mDVDY(0.0f),
        mIsUpdated(false) {}

    /* 
     * get the matrix convert world space vector to shading space
     * formed by tangent, bitangent, normal
     * normal is Vector3(0, 0, 1) in this shading space
     */
    Matrix3 Fragment::getWorldToShade() const {
        if(!mIsUpdated) {
            Vector3 n = mNormal;
            Vector3 t = normalize(mDPDU - n * dot(mDPDU, n));
            Vector3 b = cross(n, t);
            mWorldToShade = Matrix3(
                t.x, t.y, t.z,
                b.x, b.y, b.z,
                n.x, n.y, n.z);
            mIsUpdated = true;
        }
        return mWorldToShade;
    }

    void Fragment::transform(const Transform& t) {
        mPosition = t.onPoint(mPosition);
        mNormal = normalize(t.onNormal(mNormal));
        mDPDU = t.onVector(mDPDU);
        mDPDV = t.onVector(mDPDV);
        mIsUpdated = false;
    }

    size_t Geometry::nextGeometryId = 0;
    
    Geometry::Geometry(): mGeometryId(nextGeometryId++) {}

    Vector3 Geometry::sample(const Vector3& p, float u1, float u2, 
        Vector3* normal) const {
        return sample(u1, u2, normal);
    }

    float Geometry::pdf(const Vector3& p, const Vector3& wi) const {
        Ray ray(p, wi, 1e-3f);
        float epsilon;
        Fragment fragment;
        if(!intersect(ray, &epsilon, &fragment)) {
            return 0.0f;
        } 
        /* 
         * dw = da * cos(theta) / r^2
         * pdf = 1 / dw = r^2 / da * cos(theta) 
         * where r is distance between p and intersection point
         * on geometry, and theta is angle between -wi and surface normal 
         */
        float pdf = squaredLength(p - fragment.getPosition()) / 
            (area() * absdot(-wi, fragment.getNormal()));
        if(isinf(pdf)) {
            pdf = 0.0f;
        }
        return pdf;
    }


    const size_t Geometry::getId() const { return mGeometryId; }
}
