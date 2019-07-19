#ifndef GOBLIN_VECTOR_H
#define GOBLIN_VECTOR_H

#include <cassert>
#include <cmath>
#include <iostream>

namespace Goblin {
class Vector2 {
public:
    static const Vector2 UnitX;
    static const Vector2 UnitY;
    static const Vector2 Zero;
    float x, y;
public:                              
    Vector2();
    Vector2(float x, float y);
    Vector2 operator+(const Vector2& rhs) const;        
    Vector2& operator+=(const Vector2& rhs);
    Vector2 operator-(const Vector2& rhs) const;
    Vector2& operator-=(const Vector2& rhs);
    Vector2 operator*(float s) const;
    Vector2& operator*=(float s);
    Vector2 operator/(float s) const;
    Vector2& operator/=(float s);
    Vector2 operator-() const;
    const float& operator[](int i) const;
    float& operator[] (int i);
    bool operator==(const Vector2& rhs) const;
    bool operator!=(const Vector2& rhs) const;
    const float* ptr() const;
    void normalize();
    float squaredLength() const;
};

inline Vector2::Vector2() {}

inline Vector2::Vector2(float x, float y) : x(x), y(y) {}

inline Vector2 Vector2::operator+(const Vector2& rhs) const {
    return Vector2(x + rhs.x, y + rhs.y);
}

inline Vector2& Vector2::operator+=(const Vector2& rhs) {
    x += rhs.x;
    y += rhs.y;
    return *this;
}

inline Vector2 Vector2::operator-(const Vector2& rhs) const {
    return Vector2(x - rhs.x, y - rhs.y);
}

inline Vector2& Vector2::operator-=(const Vector2& rhs) {
    x -= rhs.x;
    y -= rhs.y;
    return *this;
}

inline Vector2 Vector2::operator*(float s) const {
    return Vector2(x * s, y * s);
}
        
inline Vector2& Vector2::operator*=(float s) {
    x *= s;
    y *= s;
    return *this;
}

inline Vector2 Vector2::operator/(float s) const {
    float inv = 1.0f / s;
    return Vector2(x * inv, y * inv);
}

inline Vector2& Vector2::operator/=(float s) {
    float inv = 1.0f / s;
    x *= inv;
    y *= inv;
    return *this;
}

inline Vector2 Vector2::operator-() const {
    return Vector2(-x, -y);
}
        
// to avoid const assign case like:
// const Vector2 a(1, 2);
// a[1] = 0;
inline const float& Vector2::operator[](int i) const {
    assert(i >= 0 && i < 2);
    return (&x)[i];
}

inline float& Vector2::operator[] (int i) {
    assert(i >= 0 && i < 2);
    return (&x)[i];
}

inline bool Vector2::operator==(const Vector2& rhs) const {
    return x == rhs.x && y == rhs.y;
}

inline bool Vector2::operator!=(const Vector2& rhs) const {
    return !operator==(rhs);
}

inline const float* Vector2::ptr() const {
    return &x;
}

inline float Vector2::squaredLength() const {
    return x * x + y * y;
}

inline float dot(const Vector2& lhs, const Vector2& rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y;
}

inline float absdot(const Vector2& lhs, const Vector2& rhs) {
    return fabs(dot(lhs, rhs));
}

inline float squaredLength(const Vector2& v) {
    return v.squaredLength();
}

inline float length(const Vector2& v) {
    return sqrt(squaredLength(v));
}

inline Vector2 normalize(const Vector2& v) {
    return v / length(v);
}

inline Vector2 operator*(float s, const Vector2& rhs) {
    return rhs * s;
}
    
inline std::ostream& operator<<(std::ostream& os, const Vector2& v) {
    os << "Vector2(" << v.x << ", " << v.y << ")";
    return os; 
}

class Vector3 {
public:
    static const Vector3 UnitX;
    static const Vector3 UnitY;
    static const Vector3 UnitZ;
    static const Vector3 Zero;
    float x, y, z;
public:
    Vector3();
    Vector3(float x, float y, float z);
    Vector3 operator+(const Vector3& rhs) const;
    Vector3& operator+=(const Vector3& rhs);
    Vector3 operator-(const Vector3& rhs) const;
    Vector3& operator-=(const Vector3& rhs);
    Vector3 operator*(float s) const;
    Vector3& operator*=(float s);
    Vector3 operator/(float s) const;
    Vector3& operator/=(float s);
    Vector3 operator-() const;
    const float& operator[](int i) const;
    float& operator[](int i);
    bool operator==(const Vector3& rhs) const;
    bool operator!=(const Vector3& rhs) const;
    const float* ptr() const;
    void normalize();
    float squaredLength() const;
};

inline Vector3::Vector3() {}

inline Vector3::Vector3(float x, float y, float z) : x(x), y(y), z(z) {}

inline Vector3 Vector3::operator+(const Vector3& rhs) const {
    return Vector3(x + rhs.x, y + rhs.y, z + rhs.z);
}

inline Vector3& Vector3::operator+=(const Vector3& rhs) {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    return *this;
}

inline Vector3 Vector3::operator-(const Vector3& rhs) const {
    return Vector3(x - rhs.x, y - rhs.y, z - rhs.z);
}

inline Vector3& Vector3::operator-=(const Vector3& rhs) {
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    return *this;
}

inline Vector3 Vector3::operator*(float s) const {
    return Vector3(x * s, y * s, z * s);
}

inline Vector3& Vector3::operator*=(float s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
}

inline Vector3 Vector3::operator/(float s) const {
    float inv = 1.0f / s;
    return Vector3(x * inv, y * inv, z * inv);
}
        
inline Vector3& Vector3::operator/=(float s) {
    float inv = 1.0f / s;
    x *= inv;
    y *= inv;
    z *= inv;
    return *this;
} 

inline Vector3 Vector3::operator-() const {
    return Vector3(-x, -y, -z);
}

inline const float& Vector3::operator[](int i) const {
    assert(i >= 0 && i < 3);
    return (&x)[i];
}

inline float& Vector3::operator[](int i) {
    assert(i >= 0 && i < 3);
    return (&x)[i];
}

inline bool Vector3::operator==(const Vector3& rhs) const {
    return x == rhs.x && y== rhs.y && z == rhs.z;
}

inline bool Vector3::operator!=(const Vector3& rhs) const {
    return !operator==(rhs);
}

inline const float* Vector3::ptr() const {
    return &x;
}

inline float Vector3::squaredLength() const {
    return x * x + y * y + z * z;
}

inline float dot(const Vector3& lhs, const Vector3& rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

inline float absdot(const Vector3& lhs, const Vector3& rhs) {
    return fabs(dot(lhs, rhs));
}

inline Vector3 cross(const Vector3& lhs, const Vector3& rhs) {
    return Vector3(
            lhs.y * rhs.z - lhs.z * rhs.y,
            lhs.z * rhs.x - lhs.x * rhs.z,
            lhs.x * rhs.y - lhs.y * rhs.x);
}

inline float squaredLength(const Vector3& v) {
    return v.squaredLength();
}

inline float length(const Vector3& v) {
    return sqrt(squaredLength(v));
}

inline Vector3 normalize(const Vector3& v) {
    return v / length(v);
}

inline Vector3 operator*(float s, const Vector3& rhs) {
    return rhs * s;
}

inline std::ostream& operator<<(std::ostream& os, const Vector3& v) {
    os << "Vector3(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
}

class Vector4 {
public:
    static const Vector4 UnitX;
    static const Vector4 UnitY;
    static const Vector4 UnitZ;
    static const Vector4 UnitW;
    static const Vector4 Zero;
    float x, y, z, w;
public:
    Vector4();
    Vector4(float x, float y, float z, float w);
    Vector4(const Vector3& v, float w);
    Vector4 operator+(const Vector4& rhs) const;
    Vector4& operator+=(const Vector4& rhs);
    Vector4 operator-(const Vector4& rhs) const;
    Vector4& operator-=(const Vector4& rhs);
    Vector4 operator*(float s) const;
    Vector4& operator*=(float s);
    Vector4 operator/(float s) const;
    Vector4& operator/=(float s);
    Vector4 operator-() const;
    const float& operator[](int i) const;
    float& operator[](int i);
    bool operator==(const Vector4& rhs) const;
    bool operator!=(const Vector4& rhs) const;
    const float* ptr() const;
    void normalize();
    float squaredLength() const;
};

inline Vector4::Vector4() {}

inline Vector4::Vector4(float x, float y, float z, float w) : 
    x(x), y(y), z(z), w(w) {}
        
inline Vector4::Vector4(const Vector3& v, float w) : 
    x(v.x), y(v.y), z(v.z), w(w) {}

inline Vector4 Vector4::operator+(const Vector4& rhs) const {
    return Vector4(x + rhs.x, y + rhs.y, z + rhs.z, w + rhs.w);
}

inline Vector4& Vector4::operator+=(const Vector4& rhs) {
    x += rhs.x; 
    y += rhs.y;
    z += rhs.z;
    w += rhs.w;
    return *this;
}

inline Vector4 Vector4::operator-(const Vector4& rhs) const {
    return Vector4(x - rhs.x, y - rhs.y, z - rhs.z, w - rhs.w);
}

inline Vector4& Vector4::operator-=(const Vector4& rhs) {
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    w -= rhs.w;
    return *this;
}

inline Vector4 Vector4::operator*(float s) const {
    return Vector4(x * s, y * s, z * s, w * s);
}

inline Vector4& Vector4::operator*=(float s) {
    x *= s;
    y *= s;
    z *= s;
    w *= s;
    return *this;
}

inline Vector4 Vector4::operator/(float s) const {
    float inv = 1.0f / s;
    return Vector4(x * inv, y * inv, z * inv, w * inv);
}

inline Vector4& Vector4::operator/=(float s) {
    float inv = 1.0f / s;
    x *= inv;
    y *= inv;
    z *= inv;
    w *= inv;
    return *this;
}

inline Vector4 Vector4::operator-() const {
    return Vector4(-x, -y, -z, -w);
}

inline const float& Vector4::operator[](int i) const {
    assert(i >= 0 && i < 4);
    return (&x)[i];
}

inline float& Vector4::operator[](int i) {
    assert(i >= 0 && i < 4);
    return (&x)[i];
}

inline bool Vector4::operator==(const Vector4& rhs) const {
    return x == rhs.x && y == rhs.y && z == rhs.z && w == rhs.w;
}

inline bool Vector4::operator!=(const Vector4& rhs) const {
    return !operator==(rhs);
}

inline const float* Vector4::ptr() const {
    return &x;
}

inline float Vector4::squaredLength() const {
    return x * x + y * y + z * z + w * w;
}

// return a z divide 3d vector based on the input 4d vector
inline Vector3 project(const Vector4& v) {
    float invw = v.w == 0.0f ? 1.0f : 1.0f / v.w;
    return Vector3(v.x * invw, v.y * invw, v.z * invw);
}

inline float dot(const Vector4& lhs, const Vector4& rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z + lhs.w * rhs.w;
}

inline float absdot(const Vector4& lhs, const Vector4& rhs) {
    return fabs(dot(lhs, rhs));
}

inline float squaredLength(const Vector4& v) {
    return v.squaredLength();
}

inline float length(const Vector4& v) {
    return sqrt(squaredLength(v));
}

inline Vector4 normalize(const Vector4& v) {
    return v / length(v);
}

inline Vector4 operator*(float s, const Vector4& rhs) {
    return rhs * s;
}

inline std::ostream& operator<<(std::ostream& os, const Vector4& v) {
    os << "Vector4(" << v.x << ", " << v.y << ", " << v.z << 
        ", " << v.w << ")";
    return os;
}
}
#endif //GOBLIN_VECTOR_H
