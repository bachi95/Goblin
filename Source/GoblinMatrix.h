#ifndef GOBLIN_MATRIX_H
#define GOBLIN_MATRIX_H

#include <cassert>

namespace Goblin {
    class Vector3;
    class Vector4;
    
    class Matrix3 {
    public:
        static const Matrix3 Identity;
        static const Matrix3 Zero;
    public:
        Matrix3() {}
        Matrix3(float f[9]);
        Matrix3(float m00, float m01, float m02,
                float m10, float m11, float m12,
                float m20, float m21, float m22);
        Matrix3 operator+(const Matrix3& rhs) const;
        Matrix3& operator+=(const Matrix3& rhs);
        Matrix3 operator-(const Matrix3& rhs) const;
        Matrix3& operator-=(const Matrix3& rhs);
        Matrix3 operator*(const Matrix3& rhs) const;
        Matrix3& operator*=(const Matrix3& rhs);
        Matrix3 operator*(float f) const;
        Matrix3& operator*=(float f);
        Vector3 operator*(const Vector3& rhs) const;
        Matrix3 operator/(float f) const;
        Matrix3& operator/=(float f);
        Matrix3 operator-() const; 
        
        bool operator==(const Matrix3& rhs) const;
        bool operator!=(const Matrix3& rhs) const;

        const float* operator[](int i) const {
            assert(i >= 0 && i < 3);
            return m[i];
        }

        float* operator[](int i) {
            assert(i >= 0 && i < 3);
            return m[i];
        }

        Matrix3 transpose() const;
        Matrix3& transposeSelf();

    private:
        float m[3][3];

    };

    Matrix3 operator*(float f, const Matrix3& m);
    Vector3 operator*(const Vector3& v, const Matrix3& m);
    bool inverse(Matrix3* inv, const Matrix3& m);
    bool isIdentity(const Matrix3& m);

    class Matrix4 {
    public:
        static const Matrix4 Identity;
        static const Matrix4 Zero;
    public:
        Matrix4() {}
        Matrix4(float f[16]);
        Matrix4(float m00, float m01, float m02, float m03,
                float m10, float m11, float m12, float m13,
                float m20, float m21, float m22, float m23,
                float m30, float m31, float m32, float m33);
        Matrix4 operator+(const Matrix4& rhs) const;
        Matrix4& operator+=(const Matrix4& rhs);
        Matrix4 operator-(const Matrix4& rhs) const;
        Matrix4& operator-=(const Matrix4& rhs);
        Matrix4 operator*(const Matrix4& rhs) const;
        Matrix4& operator*=(const Matrix4& rhs);
        Matrix4 operator*(float f) const;
        Matrix4& operator*=(float f);
        Vector4 operator*(const Vector4& v) const;
        Matrix4 operator/(float f) const;
        Matrix4& operator/=(float f);
        Matrix4 operator-() const;

        bool operator==(const Matrix4& rhs) const;
        bool operator!=(const Matrix4& rhs) const;

        const float* operator[](int i) const {
            assert(i >= 0 && i < 4);
            return m[i];
        }

        float* operator[](int i) {
            assert(i >= 0 && i < 4);
            return m[i];
        }

        Matrix4 transpose() const;
        Matrix4& transposeSelf();

    private:
        float m[4][4];
    };

    Matrix4 operator*(float f, const Matrix4& m);
    Vector4 operator*(const Vector4& v, const Matrix4& m);
    bool inverse(Matrix4* inv, const Matrix4& m);
    bool isIdentity(const Matrix4& m);
    
    Matrix4 matrixTranslate(const Vector3& v);
    Matrix4 matrixScale(const Vector3& s);
    //right hand coordinate, counter clockwise
    Matrix4 matrixRotationX(float angle);
    Matrix4 matrixRotationY(float angle);
    Matrix4 matrixRotationZ(float angle);
    Matrix4 matrixRotationAxis(const Vector3& v, float angle);

    Matrix4 matrixPerspectiveRHGL(float fovY, float aspect, float zn, float zf);
    Matrix4 matrixPerspectiveLHGL(float fovY, float aspect, float zn, float zf);
    Matrix4 matrixPerspectiveRHD3D(float fovY, float aspect, float zn, float zf);
    Matrix4 matrixPerspectiveLHD3D(float fovY, float aspect, float zn, float zf);

    Matrix4 matrixLookAtLH(const Vector3& eye, const Vector3& target, const Vector3& up);
    Matrix4 matrixLookAtRH(const Vector3& eye, const Vector3& target, const Vector3& up);
}

#endif //GOBLIN_MATRIX_H
