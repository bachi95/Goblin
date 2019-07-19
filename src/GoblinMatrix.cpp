#include "GoblinMatrix.h"
#include "GoblinVector.h"
#include <cstring>
#include <cmath>

namespace Goblin {
    const float MATRIX_INVERSE_EPSILON = 1e-14f;
    const float MATRIX_EPSILON = 1e-5f;
    const Matrix3 Matrix3::Identity = Matrix3(1.0f, 0.0f, 0.0f,
                                              0.0f, 1.0f, 0.0f,
                                              0.0f, 0.0f, 1.0f);

    const Matrix3 Matrix3::Zero = Matrix3(0.0f, 0.0f, 0.0f,
                                          0.0f, 0.0f, 0.0f,
                                          0.0f, 0.0f, 0.0f);

    const Matrix4 Matrix4::Identity = Matrix4(1.0f, 0.0f, 0.0f, 0.0f,
                                              0.0f, 1.0f, 0.0f, 0.0f,
                                              0.0f, 0.0f, 1.0f, 0.0f,
                                              0.0f, 0.0f, 0.0f, 1.0f);
    
    const Matrix4 Matrix4::Zero = Matrix4(0.0f, 0.0f, 0.0f, 0.0f,
                                          0.0f, 0.0f, 0.0f, 0.0f,
                                          0.0f, 0.0f, 0.0f, 0.0f,
                                          0.0f, 0.0f, 0.0f, 0.0f);
    Matrix3::Matrix3(float f[9]) {
        memcpy(m, f, sizeof(m));
    }

    Matrix3::Matrix3(float m00, float m01, float m02,
                     float m10, float m11, float m12,
                     float m20, float m21, float m22) {
        m[0][0] = m00;
        m[0][1] = m01;
        m[0][2] = m02;

        m[1][0] = m10;
        m[1][1] = m11;
        m[1][2] = m12;

        m[2][0] = m20;
        m[2][1] = m21;
        m[2][2] = m22;
    }

    Matrix3 Matrix3::operator+(const Matrix3& rhs) const {
        Matrix3 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result[i][j] = m[i][j] + rhs[i][j]; 
            }
        }
        return result;
    }

    Matrix3& Matrix3::operator+=(const Matrix3& rhs) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                m[i][j] += rhs[i][j];
            }
        }
        return *this;
    }

    Matrix3 Matrix3::operator-(const Matrix3& rhs) const {
        Matrix3 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result[i][j] = m[i][j] - rhs[i][j];
            }
        }
        return result;
    }

    Matrix3& Matrix3::operator-=(const Matrix3& rhs) {
        for (int i= 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                m[i][j] -= rhs[i][j];
            }
        }
        return *this;
    }

    Matrix3 Matrix3::operator*(const Matrix3& rhs) const {
        Matrix3 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result[i][j] = m[i][0] * rhs[0][j] + m[i][1] * rhs[1][j] +
                    m[i][2] * rhs[2][j];
            }
        }
        return result;
    }

    Matrix3& Matrix3::operator*=(const Matrix3& rhs) {
        *this = operator*(rhs);
        return *this;
    }

    Matrix3 Matrix3::operator*(float f) const {
        Matrix3 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result[i][j] = f * m[i][j];
            }
        }
        return result;
    }

    Matrix3& Matrix3::operator*=(float f) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                m[i][j] *= f;
            }
        }
        return *this;
    }
     
    Vector3 Matrix3::operator*(const Vector3& rhs) const {
        Vector3 result;
        for (int i = 0; i < 3; ++i) {
            result[i] = m[i][0] * rhs[0] + m[i][1] * rhs[1] + m[i][2] * rhs[2]; 
        }
        return result;
    }
    
    Matrix3 Matrix3::operator/(float f) const {
        float inv = 1.0f / f;
        Matrix3 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result[i][j] = m[i][j] * inv;
            }
        }
        return result;
    }

    Matrix3& Matrix3::operator/=(float f) {
        float inv = 1.0f / f;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                m[i][j] *= inv;
            }
        }
        return *this;
    }

    Matrix3 Matrix3::operator-() const {
        Matrix3 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result[i][j] = -m[i][j];
            }
        }
        return result;
    }

    bool Matrix3::operator==(const Matrix3& rhs) const {
        return memcmp(m[0], rhs[0], sizeof(m)) == 0;
    } 
    
    bool Matrix3::operator!=(const Matrix3& rhs) const {
        return !operator==(rhs);
    }

    Matrix3 Matrix3::transpose() const {
        return Matrix3(m[0][0], m[1][0], m[2][0],
                       m[0][1], m[1][1], m[2][1],
                       m[0][2], m[1][2], m[2][2]);
    }

    Matrix3& Matrix3::transposeSelf() {
        float tmp01 = m[0][1];
        m[0][1] = m[1][0];
        m[1][0] = tmp01; 
        float tmp02 = m[0][2];
        m[0][2] = m[2][0];
        m[2][0] = tmp02;
        float tmp12 = m[1][2];
        m[1][2] = m[2][1];
        m[2][1] = tmp12;
        return *this;
    }

    Matrix3 operator*(float f, const Matrix3& m) {
        return m * f;
    }

    Vector3 operator*(const Vector3& v, const Matrix3& m) {
        Vector3 result;
        for (int i = 0; i < 3; ++i) {
            result[i] = v[0] * m[0][i] + v[1] * m[1][i] + v[2] * m[2][i];
        }
        return result;
    }

    bool inverse(Matrix3* inv, const Matrix3& m) {
        assert(inv);
        (*inv)[0][0] = m[1][1] * m[2][2] - m[1][2] * m[2][1];
        (*inv)[1][0] = m[1][2] * m[2][0] - m[1][0] * m[2][2];
        (*inv)[2][0] = m[1][0] * m[2][1] - m[1][1] * m[2][0];

        float det = m[0][0] * (*inv)[0][0] + m[0][1] * (*inv)[1][0] + 
            m[0][2] * (*inv)[2][0];
        if (fabs(det) < MATRIX_INVERSE_EPSILON)
            return false;
        float invDet = 1.0f / det;

        (*inv)[0][1] = m[0][2] * m[2][1] - m[0][1] * m[2][2];
        (*inv)[0][2] = m[0][1] * m[1][2] - m[0][2] * m[1][1];
        (*inv)[1][1] = m[0][0] * m[2][2] - m[0][2] * m[2][0];
        (*inv)[1][2] = m[0][2] * m[1][0] - m[0][0] * m[1][2];
        (*inv)[2][1] = m[0][1] * m[2][0] - m[0][0] * m[2][1];
        (*inv)[2][2] = m[0][0] * m[1][1] - m[0][1] * m[1][0];

        (*inv)[0][0] *= invDet;
        (*inv)[0][1] *= invDet;
        (*inv)[0][2] *= invDet;
        (*inv)[1][0] *= invDet;
        (*inv)[1][1] *= invDet;
        (*inv)[1][2] *= invDet;
        (*inv)[2][0] *= invDet;
        (*inv)[2][1] *= invDet;
        (*inv)[2][2] *= invDet;
        
        return true;
    }

    bool isIdentity(const Matrix3& m) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (fabs(m[i][j] - float(i == j)) > MATRIX_EPSILON)
                   return false; 
            }
        }
        return true;
    }

    Matrix4::Matrix4(float f[16]) {
        memcpy(m, f, sizeof(m));
    }

    Matrix4::Matrix4(float m00, float m01, float m02, float m03,
                     float m10, float m11, float m12, float m13,
                     float m20, float m21, float m22, float m23,
                     float m30, float m31, float m32, float m33) {
        m[0][0] = m00;
        m[0][1] = m01;
        m[0][2] = m02;
        m[0][3] = m03;

        m[1][0] = m10;
        m[1][1] = m11;
        m[1][2] = m12;
        m[1][3] = m13;

        m[2][0] = m20;
        m[2][1] = m21;
        m[2][2] = m22;
        m[2][3] = m23;

        m[3][0] = m30;
        m[3][1] = m31;
        m[3][2] = m32;
        m[3][3] = m33;
    }

    Matrix4 Matrix4::operator+(const Matrix4& rhs) const {
        Matrix4 result;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                result[i][j] = m[i][j] + rhs[i][j];
            }
        }
        return result;
    }

    Matrix4& Matrix4::operator+=(const Matrix4& rhs) {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                m[i][j] += rhs[i][j];
            }
        }
        return *this;
    }

    Matrix4 Matrix4::operator-(const Matrix4& rhs) const {
        Matrix4 result;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                result[i][j] = m[i][j] - rhs[i][j];
            }
        }
        return result;
    }

    Matrix4& Matrix4::operator-=(const Matrix4& rhs) {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                m[i][j] -= rhs[i][j];
            }
        }
        return *this;
    }

    Matrix4 Matrix4::operator*(const Matrix4& rhs) const {
        Matrix4 result;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                result[i][j] = m[i][0] * rhs[0][j] + m[i][1] * rhs[1][j] +
                    m[i][2] * rhs[2][j] + m[i][3] * rhs[3][j];
            }
        }
        return result;
    }

    Matrix4& Matrix4::operator*=(const Matrix4& rhs) {
        *this = operator*(rhs);
        return *this;
    }

    Matrix4 Matrix4::operator*(float f) const {
        Matrix4 result;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                result[i][j] = f * m[i][j];
            }
        }
        return result;
    }

    Matrix4& Matrix4::operator*=(float f) {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                m[i][j] *= f;
            }
        }
        return *this;
    }

    Vector4 Matrix4::operator*(const Vector4& rhs) const {
        Vector4 result;
        for (int i = 0; i < 4; ++i) {
           result[i] = m[i][0] * rhs[0] + m[i][1] * rhs[1] + m[i][2] * rhs[2] +
               m[i][3] * rhs[3];
        }
        return result;
    }

    Matrix4 Matrix4::operator/(float f) const {
        Matrix4 result;
        float inv = 1.0f / f;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                result[i][j] = m[i][j] * inv;
            }
        }
        return result; 
    }

    Matrix4& Matrix4::operator/=(float f) {
        float inv = 1.0f / f;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                m[i][j] *= inv;
            }
        }
        return *this;
    }

    Matrix4 Matrix4::operator-() const {
        Matrix4 result;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                result[i][j] = -m[i][j];
            }
        }
        return result;
    }

    Matrix4 Matrix4::transpose() const {
        return Matrix4(m[0][0], m[1][0], m[2][0], m[3][0],
                       m[0][1], m[1][1], m[2][1], m[3][1],
                       m[0][2], m[1][2], m[2][2], m[3][2],
                       m[0][3], m[1][3], m[2][3], m[3][3]);
    }

    Matrix4& Matrix4::transposeSelf() {
        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) {
                float tmp = m[i][j];
                m[i][j] = m[j][i];
                m[j][i] = tmp;
            }
        }
        return *this;
    }
    
    bool Matrix4::operator==(const Matrix4& rhs) const {
        return memcmp(m[0], rhs[0], sizeof(m)) == 0;
    }

    bool Matrix4::operator!=(const Matrix4& rhs) const {
        return !operator==(rhs);
    }

    Matrix4 operator*(float f, const Matrix4& m) {
        return m * f;
    }

    Vector4 operator*(const Vector4& v, const Matrix4& m) {
        Vector4 result;
        for (int i = 0; i < 4; ++i) {
            result[i] = v[0] * m[0][i] + v[1] * m[1][i] + v[2] * m[2][i] +
                v[3] * m[3][i];
        }
        return result;
    }

    bool inverse(Matrix4* inv, const Matrix4& m) {
        assert(inv);
        float m00 = m[0][0];
        float m01 = m[0][1];
        float m02 = m[0][2];
        float m03 = m[0][3];
        float m10 = m[1][0];
        float m11 = m[1][1];
        float m12 = m[1][2];
        float m13 = m[1][3];
        float m20 = m[2][0];
        float m21 = m[2][1];
        float m22 = m[2][2];
        float m23 = m[2][3];
        float m30 = m[3][0];
        float m31 = m[3][1];
        float m32 = m[3][2];
        float m33 = m[3][3];

        float R23C23 = m22 * m33 - m23 * m32;
        float R23C13 = m21 * m33 - m23 * m31;
        float R23C12 = m21 * m32 - m22 * m31;
        float R23C03 = m20 * m33 - m23 * m30;
        float R23C02 = m20 * m32 - m22 * m30;
        float R23C01 = m20 * m31 - m21 * m30;

        (*inv)[0][0] = +(m11 * R23C23 - m12 * R23C13 + m13 * R23C12);
        (*inv)[1][0] = -(m10 * R23C23 - m12 * R23C03 + m13 * R23C02);
        (*inv)[2][0] = +(m10 * R23C13 - m11 * R23C03 + m13 * R23C01);
        (*inv)[3][0] = -(m10 * R23C12 - m11 * R23C02 + m12 * R23C01);

        float det = m00 * (*inv)[0][0] + m01 * (*inv)[1][0] + 
                    m02 * (*inv)[2][0] + m03 * (*inv)[3][0];
        if (fabs(det) < MATRIX_EPSILON)
            return false;
        float invDet = 1.0f / det;

        (*inv)[0][1] = -(m01 * R23C23 - m02 * R23C13 + m03 * R23C12);
        (*inv)[1][1] = +(m00 * R23C23 - m02 * R23C03 + m03 * R23C02);
        (*inv)[2][1] = -(m00 * R23C13 - m01 * R23C03 + m03 * R23C01);
        (*inv)[3][1] = +(m00 * R23C12 - m01 * R23C02 + m02 * R23C01);

        float R13C23 = m12 * m33 - m13 * m32;
        float R13C13 = m11 * m33 - m13 * m31;
        float R13C12 = m11 * m32 - m12 * m31;
        float R13C03 = m10 * m33 - m13 * m30;
        float R13C02 = m10 * m32 - m12 * m30;
        float R13C01 = m10 * m31 - m11 * m30;

        (*inv)[0][2] = +(m01 * R13C23 - m02 * R13C13 + m03 * R13C12);
        (*inv)[1][2] = -(m00 * R13C23 - m02 * R13C03 + m03 * R13C02);
        (*inv)[2][2] = +(m00 * R13C13 - m01 * R13C03 + m03 * R13C01);
        (*inv)[3][2] = -(m00 * R13C12 - m01 * R13C02 + m02 * R13C01);

        float R12C23 = m12 * m23 - m13 * m22;
        float R12C13 = m11 * m23 - m13 * m21;
        float R12C12 = m11 * m22 - m12 * m21;
        float R12C03 = m10 * m23 - m13 * m20;
        float R12C02 = m10 * m22 - m12 * m20;
        float R12C01 = m10 * m21 - m11 * m20;

        (*inv)[0][3] = -(m01 * R12C23 - m02 * R12C13 + m03 * R12C12);
        (*inv)[1][3] = +(m00 * R12C23 - m02 * R12C03 + m03 * R12C02);
        (*inv)[2][3] = -(m00 * R12C13 - m01 * R12C03 + m03 * R12C01); 
        (*inv)[3][3] = +(m00 * R12C12 - m01 * R12C02 + m02 * R12C01);

        (*inv) *= invDet;
        return true;
    }

    bool isIdentity(const Matrix4& m) {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (fabs(m[i][j] - float(i == j)) > MATRIX_EPSILON)
                   return false; 
            }
        }
        return true;
    }

    Matrix4 matrixTranslate(const Vector3& v) {
        return Matrix4(1.0f, 0.0f, 0.0f, v.x,
                       0.0f, 1.0f, 0.0f, v.y,
                       0.0f, 0.0f, 1.0f, v.z,
                       0.0f, 0.0f, 0.0f, 1.0f);
    }

    Matrix4 matrixScale(const Vector3& s) {
        return Matrix4( s.x, 0.0f, 0.0f, 0.0f,
                       0.0f,  s.y, 0.0f, 0.0f,
                       0.0f, 0.0f,  s.z, 0.0f,
                       0.0f, 0.0f, 0.0f, 1.0f);
    } 

    Matrix4 matrixRotationX(float angle) {
        float sinT = sin(angle); 
        float cosT = cos(angle);
        return Matrix4(1.0f,  0.0f,  0.0f, 0.0f,
                       0.0f,  cosT, -sinT, 0.0f,
                       0.0f,  sinT,  cosT, 0.0f,
                       0.0f,  0.0f,  0.0f, 1.0f);
    }

    Matrix4 matrixRotationY(float angle) {
        float sinT = sin(angle);
        float cosT = cos(angle);
        return Matrix4( cosT, 0.0f, sinT, 0.0f,
                        0.0f, 1.0f, 0.0f, 0.0f,
                       -sinT, 0.0f, cosT, 0.0f,
                        0.0f, 0.0f, 0.0f, 1.0f);
    }

    Matrix4 matrixRotationZ(float angle) {
        float sinT = sin(angle);
        float cosT = cos(angle);
        return Matrix4(cosT, -sinT, 0.0f, 0.0f,
                       sinT,  cosT, 0.0f, 0.0f,
                       0.0f,  0.0f, 1.0f, 0.0f,
                       0.0f,  0.0f, 0.0f, 1.0f);
    }

    // suppose we gonna rotate P around A angle theta
    // 1. decompose OP to projP(parallel to A, not affect rotation)
    // and perP(perpendicular to A, affect rotation)
    // projP = dot(A, P) * A perpP = P - dot(A, P) * A
    // 2. compose the result by linear combine perP and cross(A, P)
    // (the later one has the same length of perP and orthognal to perP and A)
    // P' = perP * cosT + cross(A, P) * sinT ->
    // P' = P * cosT + cross(A, P) * sinT + A * dot(A, P) * (1 - cosT)
    // 3, convert the above transform to matrix form:
    // R = [  c + (1 - c)AxAx  (1 - c)AxAy - sAz  (1 - c)AxAz + sAy]
    //     [(1 - c)AxAy + sAz    c + (1 - c)AyAy  (1 - c)AyAz - sAx]
    //     [(1 - c)AxAz - sAy  (1 - c)AyAz + sAx    c + (1 - c)AzAz] 
    Matrix4 matrixRotationAxis(const Vector3& v, float angle) {
        Vector3 A = normalize(v);
        float Ax = A.x;
        float Ay = A.y;
        float Az = A.z;
        float s = sin(angle);
        float c = cos(angle);
        return Matrix4(     c + (1.0f - c) * Ax * Ax, 
                       (1.0f - c) * Ax * Ay - s * Az,
                       (1.0f - c) * Ax * Az + s * Ay,
                                                0.0f,
                       (1.0f - c) * Ax * Ay + s * Az,
                            c + (1.0f - c) * Ay * Ay,
                       (1.0f - c) * Ay * Az - s * Ax,
                                                0.0f,
                       (1.0f - c) * Ax * Az - s * Ay,
                       (1.0f - c) * Ay * Az + s * Ax,
                            c + (1.0f - c) * Az * Az,
                                                0.0f,
                                                0.0f,
                                                0.0f,
                                                0.0f,
                                                1.0f);

    }

    // (x, y) is the near plane projection of point P in view space:
    // x = -n / Pz * Px  y = -n / Pz * Py
    // (x', y') is the NDC space projection result:
    // x' = (x - l) * 2 / (r - l) - 1
    // y' = (y - b) * 2 / (t - b) - 1
    // x' = 2n / (r - l) * (-Px / Pz) - (r + l) / (r - l)
    // y' = 2n / (t - b) * (-Py / Pz) - (t + b) / (t - b)
    // z' = A / z + B => -1 = A / -n + B  1 = A / f + B
    // A = 2nf / (f -n)  B = (f + n) / (f - n) 
    // z' = 2nf / (f - n) * 1 / Pz + (f + n) / (f - n)
    //
    // -x'Pz = 2n / (r - l) * Px + (r + l) / (r - l) * Pz
    // -y'Pz = 2n / (t - b) * Py + (t + b) / (t - b) * Pz
    // -z'Pz = -(f + n) / (f - n) * Pz - 2nf / (f - n)
    // the rest of them are similar( RH D3D, LH GL, LH D3D)
    Matrix4 matrixPerspectiveRHGL(float fovY, float aspect, float zn, float zf) {
        float yScale = 1.0f / tan(fovY / 2.0f);
        float xScale = yScale / aspect;
        
        Matrix4 result = Matrix4::Zero;
        result[0][0] = xScale;
        result[1][1] = yScale;
        result[2][2] = (zf + zn) / (zn - zf);
        result[2][3] = 2.0f * zn * zf / (zn - zf),
        result[3][2] = -1.0f;
        return result; 
    }

    Matrix4 matrixPerspectiveLHGL(float fovY, float aspect, float zn, float zf) {
        float yScale = 1.0f / tan(fovY / 2.0f);
        float xScale = yScale / aspect;

        Matrix4 result = Matrix4::Zero;
        result[0][0] = xScale;
        result[1][1] = yScale;
        result[2][2] = (zf + zn) / (zf - zn);
        result[2][3] = 2.0f * zn * zf / (zn - zf),
        result[3][2] = 1.0f;
        return result;
    }

    Matrix4 matrixPerspectiveRHD3D(float fovY, float aspect, float zn, float zf) {
        float yScale = 1.0f / tan(fovY / 2.0f);
        float xScale = yScale / aspect;

        Matrix4 result = Matrix4::Zero;
        result[0][0] = xScale;
        result[1][1] = yScale;
        result[2][2] = zf / (zn - zf);
        result[2][3] = zn * zf / (zn - zf),
        result[3][2] = -1.0f;
        return result;
    }

    Matrix4 matrixPerspectiveLHD3D(float fovY, float aspect, float zn, float zf) {
        float yScale = 1.0f / tan(fovY / 2.0f);
        float xScale = yScale / aspect;

        Matrix4 result = Matrix4::Zero;
        result[0][0] = xScale;
        result[1][1] = yScale;
        result[2][2] = zf / (zf - zn);
        result[2][3] = -zn * zf / (zf - zn),
        result[3][2] = 1.0f;
        return result;
    }

    Matrix4 matrixOrthoRHGL(float w, float h, float zn, float zf) {
        Matrix4 result = Matrix4::Zero;
        result[0][0] = 2.0f / w;
        result[1][1] = 2.0f / h;
        result[2][2] = 2.0f / (zn - zf);
        result[2][3] = (zf + zn) / (zn -zf);
        result[3][3] = 1.0f;
        return result;
    }

    Matrix4 matrixOrthoLHGL(float w, float h, float zn, float zf) {
        Matrix4 result = Matrix4::Zero;
        result[0][0] = 2.0f / w;
        result[1][1] = 2.0f / h;
        result[2][2] = 2.0f / (zf - zn);
        result[2][3] = (zf + zn) / (zn -zf);
        result[3][3] = 1.0f;
        return result;
    }

    Matrix4 matrixOrthoRHD3D(float w, float h, float zn, float zf) {
        Matrix4 result = Matrix4::Zero;
        result[0][0] = 2.0f / w;
        result[1][1] = 2.0f / h;
        result[2][2] = 1.0f / (zn - zf);
        result[2][3] = zn / (zn -zf);
        result[3][3] = 1.0f;
        return result;
    }

    Matrix4 matrixOrthoLHD3D(float w, float h, float zn, float zf) {
        Matrix4 result = Matrix4::Zero;
        result[0][0] = 2.0f / w;
        result[1][1] = 2.0f / h;
        result[2][2] = 1.0f / (zf - zn);
        result[2][3] = zn / (zn -zf);
        result[3][3] = 1.0f;
        return result;
    }

    Matrix4 matrixLookAtLH(const Vector3& eye, const Vector3& target, const Vector3& up) {
        Vector3 zAxis = normalize(target - eye);
        Vector3 xAxis = normalize(cross(up, zAxis));
        Vector3 yAxis = cross(zAxis, xAxis);

        return Matrix4(xAxis.x, xAxis.y, xAxis.z, -dot(xAxis, eye),
                       yAxis.x, yAxis.y, yAxis.z, -dot(yAxis, eye),
                       zAxis.x, zAxis.y, zAxis.z, -dot(zAxis, eye),
                          0.0f,    0.0f,    0.0f,            1.0f);
    }

    Matrix4 matrixLookAtRH(const Vector3& eye, const Vector3& target, const Vector3& up) {
        Vector3 zAxis = normalize(eye - target);
        Vector3 xAxis = normalize(cross(up, zAxis));
        Vector3 yAxis = cross(zAxis, xAxis);

        return Matrix4(xAxis.x, xAxis.y, xAxis.z, -dot(xAxis, eye),
                       yAxis.x, yAxis.y, yAxis.z, -dot(yAxis, eye),
                       zAxis.x, zAxis.y, zAxis.z, -dot(zAxis, eye),
                          0.0f,    0.0f,    0.0f,            1.0f);
    }

}

