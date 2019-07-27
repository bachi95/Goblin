#include "GoblinQuaternion.h"
#include "GoblinMatrix.h"
#include "GoblinUtils.h"
#include <cmath>
namespace Goblin {

const Quaternion Quaternion::Identity = Quaternion(1.0f, 0.0f, 0.0f, 0.0f);

Quaternion::Quaternion(const Vector3& axis, float angle) {
    float t = angle * 0.5f;
    Vector3 unitAxis = normalize(axis);
    float sinT = sin(t);
    w = cos(t);
    v = unitAxis * sinT;
}

// this can be reverse derived from the toMatrix
// equation below, for more detail reference 
// Ken Shoemake's quaternion tutorial
// http://www.cs.ucr.edu/~vbz/resources/quatut.pdf
Quaternion::Quaternion(const Matrix3& R) {
    float q[4];
    float trace = R[0][0] + R[1][1] + R[2][2];
    if (trace > 0.0f) {
        float s = sqrt(trace + 1.0f);
        q[3] = s * 0.5f;
        float t = 0.5f / s;
        q[0] = (R[2][1] - R[1][2]) * t;
        q[1] = (R[0][2] - R[2][0]) * t;
        q[2] = (R[1][0] - R[0][1]) * t;
    } else {
        int i = 0;
        if (R[1][1] > R[0][0]) {
            i = 1;
        }
        if (R[2][2] > R[i][i]) {
            i = 2;
        }
        static const int next[3] = {1, 2, 0};
        int j = next[i];
        int k = next[j];

        float s = sqrt(R[i][i] - R[j][j] - R[k][k] + 1.0f);
        q[i] = s * 0.5f;
        float t;
        t = s != 0.0f ? 0.5f / s : s;
        q[3] = (R[k][j] - R[j][k]) * t;
        q[j] = (R[j][i] + R[i][j]) * t;
        q[k] = (R[k][i] + R[i][k]) * t;
    }
    v = Vector3(q[0], q[1], q[2]);
    w = q[3];
}

Matrix4 Quaternion::toMatrix() const {
    float x2 = 2.0f * v.x;
    float y2 = 2.0f * v.y;
    float z2 = 2.0f * v.z;
    float xx2 = x2 * v.x;
    float xy2 = x2 * v.y;
    float xz2 = x2 * v.z;
    float xw2 = x2 * w;
    float yy2 = y2 * v.y;
    float yz2 = y2 * v.z;
    float yw2 = y2 * w;
    float zz2 = z2 * v.z;
    float zw2 = z2 * w;

    return Matrix4(1 - yy2 - zz2,    xy2 - zw2,    xz2 + yw2, 0.0f,
                        xy2 + zw2, 1 -xx2 - zz2,    yz2 - xw2, 0.0f,
                        xz2 - yw2,    yz2 + xw2, 1 -xx2 - yy2, 0.0f,
                            0.0f,         0.0f,         0.0f, 1.0f);
}

// the NVIDIA optimized version of qpq^-1
// quaternion w + v => v = tA
// where A is unit vector(rotation axis)
// t = sin(theta/2) s = cos(theta/2)
// point to rotate: p
// p + 2w(cross(v, p)) + 2(cross(v, cross(v, p))) =
// p + 2wt(cross(A, p)) + 2(tA(dot(tA, p)) - (t^2)(A^2)p) =
// ( by just unroll the cross(v, cross(v, p)) )
// (1-2t^2)p + 2wt((cross(A, p)) + 2tA(dot(tA, p)) =
// (w^2 - t^2)p + 2wt(cross(A, p)) + 2tA(dot(tA, P)) =
// qpq^-1
Vector3 Quaternion::operator*(const Vector3& p) const {
    Vector3 uv = cross(v, p);
    Vector3 uuv = cross(v, uv);
    uv *= (2.0f * w);
    uuv *= 2.0f;
    return p + uv + uuv;
}

Quaternion normalize(const Quaternion& q) {
    float inv = 1.0f / sqrt(q.norm());
    Quaternion result(q);
    result.v *= inv;
    result.w *= inv;
    return result;
}

Quaternion eulerToQuaternion(const Vector3& xyzAngle,
    const std::string& rotationOrder) {
    RotationOrder order = XYZ;
    if (rotationOrder == "xyz") {
        order = XYZ;
    } else if (rotationOrder == "xzy") {
        order = XZY;
    } else if (rotationOrder == "yxz") {
        order = YXZ;
    } else if (rotationOrder == "yzx") {
        order = YZX;
    } else if (rotationOrder == "zxy") {
        order = ZXY;
    } else if (rotationOrder == "zyx") {
        order = ZYX;
    } else {
		std::cerr << "unrecognized rotation order " << rotationOrder <<
            ", fall back to XYZ" << std::endl;
    }
    return eulerToQuaternion(xyzAngle, order);
}

Quaternion eulerToQuaternion(const Vector3& xyzAngle, RotationOrder order) {
    Quaternion qx(Vector3::UnitX, radians(xyzAngle.x));
    Quaternion qy(Vector3::UnitY, radians(xyzAngle.y));
    Quaternion qz(Vector3::UnitZ, radians(xyzAngle.z));
    Quaternion result;
    if (order == XYZ) {
        result = qz * qy * qx;
    } else if (order == XZY) {
        result = qy * qz * qx;
    } else if (order == YXZ) {
        result = qz * qx * qy;
    } else if (order == YZX) {
        result = qx * qz * qy;
    } else if (order == ZXY) {
        result = qy * qx * qz;
    } else if (order == ZYX) {
        result = qx * qy * qz;
    } else {
		std::cerr << "unrecognized rotation order, fall back to XYZ" << std::endl;
        result = qz * qy * qx;
    }
    return result;
}

std::ostream& operator<<(std::ostream& os, const Quaternion& q) {
    os << "Quaternion(" << q.w << ", " << 
        q.v.x << ", " << q.v.y << ", " << q.v.z <<")";
    return os;
}
}
