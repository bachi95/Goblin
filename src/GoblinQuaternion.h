#ifndef GOBLIN_QUATERNION_H
#define GOBLIN_QUATERNION_H

#include "GoblinVector.h"

#include <string>

namespace Goblin {
class Matrix3;
class Matrix4;

enum RotationOrder {
    XYZ,
    XZY,
    YXZ,
    YZX,
    ZXY,
    ZYX
};

class Quaternion {
public:
	static const Quaternion Identity;
public:
	Quaternion() = default;

	Quaternion(float w, float x, float y, float z) : w(w), v(x, y, z) {}

	Quaternion(float w, const Vector3& v) : w(w), v(v) {}

	Quaternion(const Vector3& axis, float angle);

	Quaternion(const Matrix3& R);

	Quaternion conjugate() const {
		return Quaternion(w, -v);
	}

	float norm() const {
		return w * w + v.squaredLength();
	}

	Matrix4 toMatrix() const;

	Quaternion operator*(const Quaternion& rhs) const {
		return Quaternion(w * rhs.w - dot(v, rhs.v),
			w * rhs.v + rhs.w * v + cross(v, rhs.v));
	}

	Vector3 operator*(const Vector3& p) const;

	bool operator==(const Quaternion& rhs) const {
		return w == rhs.w && v == rhs.v;
	}

	bool operator!=(const Quaternion& rhs) const {
		return !operator==(rhs);
	}

public:
    float w;
    Vector3 v;
}; 

Quaternion normalize(const Quaternion& q);

Quaternion eulerToQuaternion(const Vector3& xyzAngle, 
    const std::string& rotationOrder);

Quaternion eulerToQuaternion(const Vector3& xyzAngle, RotationOrder order);

std::ostream& operator<<(std::ostream& os, const Quaternion& q);
}

#endif //GOBLIN_QUATERNION_H
