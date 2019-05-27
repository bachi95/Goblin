#include "GoblinTransform.h"
#include "GoblinRay.h"
#include "GoblinBBox.h"
#include <iostream>

namespace Goblin {

    Transform::Transform():
        mCachedMatrix(Matrix4::Identity), mCachedInverse(Matrix4::Identity),
        mPosition(Vector3::Zero), mOrientation(Vector3::UnitY, 0.0f),
        mScale(Vector3(1.0f, 1.0f, 1.0f)), mIsUpdated(true) {}

        Transform::Transform(const Vector3& position, 
            const Quaternion& orientation,
            const Vector3& scale): mPosition(position),
            mOrientation(orientation), mScale(scale),
            mIsUpdated(false) {
            update();
        }

        void Transform::setPosition(const Vector3& position) {
            mPosition = position;
            mIsUpdated = false;
        }

        void Transform::setOrientation(const Quaternion& orientation) {
            mOrientation = orientation;
            mIsUpdated = false;
        }

        void Transform::setScale(const Vector3& scale) {
            mScale = scale;
            mIsUpdated = false;
        }

        const Vector3& Transform::getPosition() const {
            return mPosition;
        }

        const Quaternion& Transform::getOrientation() const {
            return mOrientation;
        }

        const Matrix4& Transform::getMatrix() const {
            if (!isUpdated()) {
                update();
            }
            return mCachedMatrix;
        }

        const Matrix4& Transform::getInverse() const {
            if (!isUpdated()) {
                update();
            }
            return mCachedInverse;
        }

        const Vector3& Transform::getScale() const {
            return mScale;
        }

        void Transform::roll(float angle) {
            rotate(mOrientation * Vector3::UnitZ, angle);
        }

        void Transform::pitch(float angle) {
            rotate(mOrientation * Vector3::UnitX, angle);
        }

        void Transform::yaw(float angle) {
            rotate(mOrientation * Vector3::UnitY, angle);
        }

        void Transform::rotateX(float angle) {
            rotate(Vector3::UnitX, angle);
        }

        void Transform::rotateY(float angle) {
            rotate(Vector3::UnitY, angle);
        }

        void Transform::rotateZ(float angle) {
            rotate(Vector3::UnitZ, angle);
        }

        void Transform::rotate(const Vector3& axis, float angle) {
            mOrientation = normalize(Quaternion(axis, angle) * mOrientation);
            mIsUpdated = false;
        }

        void Transform::translate(const Vector3& d) {
            mPosition += d;
            mIsUpdated = false;
        }

        Vector3 Transform::onPoint(const Vector3& p) const {
            const Matrix4& M = getMatrix();
            return Vector3(
                M[0][0] * p.x + M[0][1] * p.y + M[0][2] * p.z + M[0][3],
                M[1][0] * p.x + M[1][1] * p.y + M[1][2] * p.z + M[1][3],
                M[2][0] * p.x + M[2][1] * p.y + M[2][2] * p.z + M[2][3]);
        }

        Vector3 Transform::onNormal(const Vector3& n) const {
            const Matrix4& invT = getInverse().transpose();
            return Vector3(
                invT[0][0] * n.x + invT[0][1] * n.y + invT[0][2] * n.z,
                invT[1][0] * n.x + invT[1][1] * n.y + invT[1][2] * n.z,
                invT[2][0] * n.x + invT[2][1] * n.y + invT[2][2] * n.z);
        }

        Vector3 Transform::onVector(const Vector3& v) const {
            const Matrix4& M = getMatrix();
            return Vector3(
                M[0][0] * v.x + M[0][1] * v.y + M[0][2] * v.z,
                M[1][0] * v.x + M[1][1] * v.y + M[1][2] * v.z,
                M[2][0] * v.x + M[2][1] * v.y + M[2][2] * v.z);
        }

        Ray Transform::onRay(const Ray& r) const {
            return Ray(onPoint(r.o), onVector(r.d), r.mint, r.maxt, r.depth);
        }

        BBox Transform::onBBox(const BBox& b) const {
            BBox rv(onPoint(b.pMin));
            rv.expand(onPoint(Vector3(b.pMax.x, b.pMin.y, b.pMin.z)));
            rv.expand(onPoint(Vector3(b.pMin.x, b.pMax.y, b.pMin.z)));
            rv.expand(onPoint(Vector3(b.pMin.x, b.pMin.y, b.pMax.z)));
            rv.expand(onPoint(Vector3(b.pMax.x, b.pMax.y, b.pMin.z)));
            rv.expand(onPoint(Vector3(b.pMax.x, b.pMin.y, b.pMax.z)));
            rv.expand(onPoint(Vector3(b.pMin.x, b.pMax.y, b.pMax.z)));
            rv.expand(onPoint(b.pMax));
            return rv;
        }

        Vector3 Transform::invertPoint(const Vector3& p) const {
            const Matrix4& M = getInverse();
            return Vector3(
                M[0][0] * p.x + M[0][1] * p.y + M[0][2] * p.z + M[0][3],
                M[1][0] * p.x + M[1][1] * p.y + M[1][2] * p.z + M[1][3],
                M[2][0] * p.x + M[2][1] * p.y + M[2][2] * p.z + M[2][3]);
        }

        Vector3 Transform::invertNormal(const Vector3& n) const {
            const Matrix4& invT = getMatrix().transpose();
            return Vector3(
                invT[0][0] * n.x + invT[0][1] * n.y + invT[0][2] * n.z,
                invT[1][0] * n.x + invT[1][1] * n.y + invT[1][2] * n.z,
                invT[2][0] * n.x + invT[2][1] * n.y + invT[2][2] * n.z);
        }

        Vector3 Transform::invertVector(const Vector3& v) const {
            const Matrix4& M = getInverse();
            return Vector3(
                M[0][0] * v.x + M[0][1] * v.y + M[0][2] * v.z,
                M[1][0] * v.x + M[1][1] * v.y + M[1][2] * v.z,
                M[2][0] * v.x + M[2][1] * v.y + M[2][2] * v.z);
        }

        Ray Transform::invertRay(const Ray& r) const {
            return Ray(invertPoint(r.o), invertVector(r.d), 
                r.mint, r.maxt, r.depth);
        }

        BBox Transform::invertBBox(const BBox& b) const {
            BBox rv(invertPoint(b.pMin));
            rv.expand(invertPoint(Vector3(b.pMax.x, b.pMin.y, b.pMin.z)));
            rv.expand(invertPoint(Vector3(b.pMin.x, b.pMax.y, b.pMin.z)));
            rv.expand(invertPoint(Vector3(b.pMin.x, b.pMin.y, b.pMax.z)));
            rv.expand(invertPoint(Vector3(b.pMax.x, b.pMax.y, b.pMin.z)));
            rv.expand(invertPoint(Vector3(b.pMax.x, b.pMin.y, b.pMax.z)));
            rv.expand(invertPoint(Vector3(b.pMin.x, b.pMax.y, b.pMax.z)));
            rv.expand(invertPoint(b.pMax));
            return rv;
        }

        bool Transform::isUpdated() const {
            return mIsUpdated;
        }

        void Transform::update() const {
            Matrix4 S = matrixScale(mScale);
            Matrix4 R = mOrientation.toMatrix();

            mCachedMatrix = R * S;
            mCachedMatrix[0][3] = mPosition.x;
            mCachedMatrix[1][3] = mPosition.y;
            mCachedMatrix[2][3] = mPosition.z;
            
            inverse(&mCachedInverse, mCachedMatrix);
            mIsUpdated = true;
        }
}