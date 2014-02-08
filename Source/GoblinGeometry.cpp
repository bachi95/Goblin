#include "GoblinGeometry.h"
#include "GoblinRay.h"

namespace Goblin {
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
        float pdf = squaredLength(p - fragment.position) / 
            (area() * absdot(-wi, fragment.normal));
        return pdf;    
    }


    const size_t Geometry::getId() const { return mGeometryId; }
}
