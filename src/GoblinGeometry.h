#ifndef GOBLIN_GEOMETRY_H
#define GOBLIN_GEOMETRY_H
#include "GoblinVertex.h"
#include "GoblinBBox.h"
#include "GoblinUtils.h"
#include "GoblinMatrix.h"

namespace Goblin {
class BBox;
class Ray;
class Transform;

class Fragment {
public:
	Fragment() :
		mDPDX(Vector3::Zero), mDPDY(Vector3::Zero),
		mDUDX(0.0f), mDVDX(0.0f), mDUDY(0.0f), mDVDY(0.0f),
		mIsUpdated(false)
	{
	}

	Fragment(const Vector3& p, const Vector3& n, const Vector2& uv,
		const Vector3& dpdu, const Vector3& dpdv);

	const Vector3& getPosition() const {
		return mPosition;
	}

	const Vector3& getNormal() const {
		return mNormal;
	}

	const Vector2& getUV() const {
		return mUV;
	}

	const Vector3& getDPDU() const {
		return mDPDU;
	}

	const Vector3& getDPDV() const {
		return mDPDV;
	}

	const Vector3& getDPDX() const {
		return mDPDX;
	}

	const Vector3& getDPDY() const {
		return mDPDY;
	}

	float getDUDX() const {
		return mDUDX;
	}

	float getDVDX() const {
		return mDVDX;
	}

	float getDUDY() const {
		return mDUDY;
	}

	float getDVDY() const {
		return mDVDY;
	}

    Matrix3 getWorldToShade() const;

    void setPosition(const Vector3& position) {
		mIsUpdated = false;
		mPosition = position;
	}

    void setNormal(const Vector3& normal) {
		mIsUpdated = false;
		mNormal = normal;
	}

    void setUV(const Vector2& uv) {
		mIsUpdated = false;
		mUV = uv;
	}

    void setDPDU(const Vector3& dpdu) {
		mIsUpdated = false;
		mDPDU = dpdu;
	}

    void setDPDV(const Vector3& dpdv) {
		mIsUpdated = false;
		mDPDV = dpdv;
	}

    void setDPDX(const Vector3& dpdx) {
		mIsUpdated = false;
		mDPDX = dpdx;
	}

    void setDPDY(const Vector3& dpdy) {
		mIsUpdated = false;
		mDPDY = dpdy;
	}

    void setUVDifferential(float dudx, float dvdx, float dudy, float dvdy) {
		mIsUpdated = false;
		mDUDX = dudx;
		mDVDX = dvdx;
		mDUDY = dudy;
		mDVDY = dvdy;
	}

    void transform(const Transform& t);

private:
    Vector3 mPosition;
    Vector3 mNormal;
    Vector2 mUV;
    Vector3 mDPDU;
    Vector3 mDPDV;
    Vector3 mDPDX;
    Vector3 mDPDY;
    float mDUDX;
    float mDVDX;
    float mDUDY;
    float mDVDY;
    mutable bool mIsUpdated;
    mutable Matrix3 mWorldToShade;
};

struct TriangleIndex {
    unsigned int v[3];
};

class Geometry;
typedef std::vector<const Geometry*> GeometryList;

class Geometry {
public:
    Geometry();

	virtual ~Geometry() = default;

	virtual bool intersectable() const { return true; }

    virtual bool intersect(const Ray& ray, float* epsilon,
        Fragment* fragment) const = 0;

	virtual bool occluded(const Ray& ray) const = 0;

    virtual Vector3 sample(float u1, float u2, Vector3* normal) const;

    virtual Vector3 sample(const Vector3& p, float u1, float u2,
        Vector3* normal) const;

    // pdf with respect to solid angle from p and extanded by this geometry
    virtual float pdf(const Vector3& p, const Vector3& wi) const;

    virtual float area() const = 0;

    virtual BBox getObjectBound() const = 0;

    virtual void refine(GeometryList& refinedGeometries) const;

    const size_t getId() const;

    static void clearGeometryCache();

protected:
    static size_t nextGeometryId;
    static std::map<size_t, Geometry*> geometryCache;
    size_t mGeometryId;
};

inline Vector3 Geometry::sample(float u1, float u2,
    Vector3* normal) const {
    std::cerr << "unimplemented Geometry::sample(float, float, Vector3*)"
        << std::endl;
    throw std::exception();
}

inline void Geometry::refine(GeometryList& refinedGeometries) const {
    std::cerr << "unimplemented Geometry::refine" << std::endl;
    throw std::exception();
}

inline void Geometry::clearGeometryCache() {
    std::map<size_t, Geometry*>::iterator it;
    for (it = geometryCache.begin(); it != geometryCache.end(); ++it) {
        delete it->second;
    }
    geometryCache.clear();
}

}

#endif //GOBLIN_GEOMETRY_H