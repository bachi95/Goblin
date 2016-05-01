#include "GoblinCamera.h"
#include "GoblinFilm.h"
#include "GoblinRay.h"
#include "GoblinSampler.h"
#include "GoblinUtils.h"

namespace Goblin{

    const Vector3 Camera::sInvalidPixel =
        Vector3(-INFINITY, -INFINITY, -INFINITY);

    Camera::Camera(const Vector3& position, const Quaternion& orientation,
        float zn, float zf, Film* film):
        mPosition(position), mOrientation(orientation),
        mZNear(zn), mZFar(zf), mFilm(film), mFilmArea(0.0f),
        mPixelArea(0.0f), mIsUpdated(false) {
        float xRes = static_cast<float>(film->getXResolution());
        float yRes = static_cast<float>(film->getYResolution());
        mAspectRatio = xRes / yRes; 
    }

    Camera::~Camera() {
        if(mFilm != NULL) {
            delete mFilm;
            mFilm = NULL;
        }
    }

    const Matrix4& Camera::getWorldMatrix() {
        if(!isUpdated()) {
            update();
        }
        return mWorld;
    }

    const Matrix4& Camera::getViewMatrix() {
        if(!isUpdated()) {
            update();
        }
        return mView;
    }

    Vector3 Camera::worldToScreen(const Vector3& pWorld) const {
        Vector4 pView = mView * Vector4(pWorld, 1.0f); 
        Vector4 pNDC = mProj * pView;
        pNDC /= pNDC.w;
        float screenX = (pNDC.x + 1.0f) * 0.5f * mFilm->getXResolution();
        float screenY = (1.0f - pNDC.y) * 0.5f * mFilm->getYResolution();
        return Vector3(screenX, screenY, pView.z);
    }

    Vector3 Camera::worldToScreen(const Vector3& pWorld,
        const Vector3& pLens) const {
        return worldToScreen(pWorld);
    }
    
    void Camera::update() {
        Matrix4 R = mOrientation.toMatrix();

        Vector3 right(R[0][0], R[1][0], R[2][0]);
        Vector3    up(R[0][1], R[1][1], R[2][1]);
        Vector3  look(R[0][2], R[1][2], R[2][2]);

        mWorld = Matrix4(
             right.x, up.x, look.x, mPosition.x,
             right.y, up.y, look.y, mPosition.y,
             right.z, up.z, look.z, mPosition.z,
                0.0f, 0.0f,   0.0f,        1.0f);

        float x = -dot(mPosition, right);
        float y = -dot(mPosition, up);
        float z = -dot(mPosition, look);

        mView = Matrix4(
            right.x, right.y,  right.z,    x,
               up.x,    up.y,     up.z,    y,
             look.x,  look.y,   look.z,    z,
               0.0f,    0.0f,     0.0f, 1.0f);

        mIsUpdated = true;
    }

    PerspectiveCamera::PerspectiveCamera(const Vector3& position, 
        const Quaternion& orientation, float fov, float zn, float zf, 
        float lensRadius, float focalDistance, Film* film):
        Camera(position, orientation, zn, zf, film),
        mFOV(fov), mLensRadius(lensRadius), mFocalDistance(focalDistance) {
        mProj = matrixPerspectiveLHD3D(mFOV, mAspectRatio, mZNear, mZFar);
        // pro[0][0] = 1 / (tan(fov / 2) * aspectRatio)
        // pro[1][1] = 1 / (tan(fov / 2))
        float filmWorldHeight = 2.0f * mFocalDistance * tan(0.5f * mFOV);
        float filmWorldWidth = filmWorldHeight * mAspectRatio;
        mFilmArea = filmWorldHeight * filmWorldWidth;
        mPixelArea = mFilmArea *
            mFilm->getInvXResolution() * mFilm->getInvYResolution();
        mFilm->setFilmArea(mFilmArea);
        update();
    }

    float PerspectiveCamera::generateRay(const Sample& sample, 
        RayDifferential* ray) const {
        float invXRes = mFilm->getInvXResolution();
        float invYRes = mFilm->getInvYResolution();
        float xNDC = +2.0f * sample.imageX * invXRes - 1.0f;
        float yNDC = -2.0f * sample.imageY * invYRes + 1.0f;
        float dxNDC = +2.0f * (sample.imageX + 1.0f) * invXRes - 1.0f;
        float dyNDC = -2.0f * (sample.imageY + 1.0f) * invYRes + 1.0f;
        // from NDC space to view space
        // xView = xNDC * zView * tan(fov / 2) * aspectRatio
        // yView = yNDC * zView * tan(fov / 2)
        // in projection matrix pro,
        // pro[0][0] = 1 / (tan(fov / 2) * aspectRatio)
        // pro[1][1] = 1 / (tan(fov / 2))
        float zView = 1.0f;
        float xView = xNDC / mProj[0][0];
        float yView = yNDC / mProj[1][1];
        Vector3 viewDir = Vector3(xView, yView, zView);

        float dxView = dxNDC / mProj[0][0];
        Vector3 dxViewDir(dxView, yView, zView);
        float dyView = dyNDC / mProj[1][1];
        Vector3 dyViewDir(xView, dyView, zView);

        if(mLensRadius == 0.0f) {
            // view space to world space
            ray->o = ray->dxOrigin = ray->dyOrigin = mPosition;
            ray->d = mOrientation * normalize(viewDir);
            ray->dxDir = mOrientation * normalize(dxViewDir);
            ray->dyDir = mOrientation * normalize(dyViewDir);
            
        } else {
            // perturb the ray direction for DOF effect
            float ft = mFocalDistance / viewDir.z;
            Vector3 pFocus = viewDir * ft;
            Vector3 pDxFocus = dxViewDir * ft;
            Vector3 pDyFocus = dyViewDir * ft;
            Vector2 lensSample = mLensRadius *
                uniformSampleDisk(sample.lensU1, sample.lensU2);
            Vector3 viewOrigin(lensSample.x, lensSample.y, 0.0f);
            ray->o = ray->dxOrigin = ray->dyOrigin = 
                mOrientation * viewOrigin + mPosition;
            ray->d = mOrientation * normalize(pFocus - viewOrigin);
            ray->dxDir = mOrientation * normalize(pDxFocus - viewOrigin);
            ray->dyDir = mOrientation * normalize(pDyFocus - viewOrigin);
        }
        ray->mint = 1e-3f;
        ray->maxt = INFINITY;
        ray->depth = 0;
        ray->hasDifferential = true;
        return 1.0f;
    }

    Vector3 PerspectiveCamera::samplePosition(const Sample& sample,
        Vector3* surfaceNormal, float* pdfArea) const {
        Vector3 pCamera;
        if (mLensRadius > 0.0f) {
            Vector2 lensSample = mLensRadius *
                uniformSampleDisk(sample.lensU1, sample.lensU2);
            Vector3 viewOrigin(lensSample.x, lensSample.y, 0.0f);
            pCamera = mOrientation * viewOrigin + mPosition;
        } else {
            pCamera = mPosition;
        }
        if (pdfArea) {
            *pdfArea = mLensRadius > 0.0f ?
                1.0f / (mLensRadius * mLensRadius * PI) : 1.0f;
        }
        *surfaceNormal = getLook();
        return pCamera;
    }

    Vector3 PerspectiveCamera::sampleDirection(const Sample& sample,
        const Vector3& pCamera, float* We, float* pdfW) const {
        float invXRes = mFilm->getInvXResolution();
        float invYRes = mFilm->getInvYResolution();
        float xNDC = +2.0f * sample.imageX * invXRes - 1.0f;
        float yNDC = -2.0f * sample.imageY * invYRes + 1.0f;
        // from NDC space to view space
        // xView = xNDC * zView * tan(fov / 2) * aspectRatio
        // yView = yNDC * zView * tan(fov / 2)
        // in projection matrix pro,
        // pro[0][0] = 1 / (tan(fov / 2) * aspectRatio)
        // pro[1][1] = 1 / (tan(fov / 2))
        float zView = 1.0f;
        float xView = xNDC / mProj[0][0];
        float yView = yNDC / mProj[1][1];
        Vector3 viewDir = Vector3(xView, yView, zView);
        Vector3 pFocusView = mFocalDistance * viewDir;
        Vector3 pFocusWorld = mOrientation * pFocusView + mPosition;
        Vector3 sampleDir = pFocusWorld - pCamera;
        float lensToFilmDistance2 = sampleDir.squaredLength();
        sampleDir.normalize();
        // let w be soilid angle start on a a point in lens and extend by film
        // Power = integrate(integrate(L *cosTheta) over w) over lensArea
        // the estimator for the above integration is:
        // L * cosTheta / (pdf(leansArea) * pdf(w)) =
        // L * cosTheta * lensArea * filmArea * cosTheta / (lensToFilm^2) =
        // L * filmArea * lensArea * G where G is the geometry term
        // expectationValue(L) = expectationValue(Power * We)
        // We = 1 / (filmArea * lensArea * G)
        float cosTheta = absdot(mOrientation * Vector3::UnitZ, sampleDir);
        float G = cosTheta * cosTheta / lensToFilmDistance2;
        float lensArea = mLensRadius * mLensRadius * PI;
        *We = mLensRadius > 0.0f ?
            1.0f / (mFilmArea * lensArea * G) : 1.0f / (mFilmArea * G);
        if (pdfW) {
            *pdfW = lensToFilmDistance2 / (mFilmArea * cosTheta);
        }
        return sampleDir;
    }

    float PerspectiveCamera::pdfPosition(const Vector3& p) const {
        return mLensRadius > 0.0f ?
            1.0f / mLensRadius * mLensRadius * PI : 0.0f;
    }

    float PerspectiveCamera::pdfDirection(const Vector3& p,
        const Vector3& wo) const {
        // this method doesn't check whether p is projected to film plane
        // based on the assumption that caller already cull out the
        // position that is out of camera frustum

        // pdfW = pdfA / G = 1 / filmArea * r * r / cosTheta
        // wehre cosTheta = dot(wo, n) and r = focalDistance / cosTheta
        // pdfW = focalDistance^2 / (filmArea * cosTheata^3)
        float cosTheta = dot(getLook(),wo);
        float pdfW = mFocalDistance * mFocalDistance /
            (mFilmArea * cosTheta * cosTheta * cosTheta);
        return pdfW;
    }

    float PerspectiveCamera::evalWe(const Vector3& pCamera,
        const Vector3& pWorld) const {
        if (worldToScreen(pWorld, pCamera) == Camera::sInvalidPixel) {
            return 0.0f;
        }
        Vector4 pView = mView * Vector4(pWorld, 1.0f);
        Vector4 pLensLocal = mView * Vector4(pCamera, 1.0f);
        Vector4 dir(pView - pLensLocal);
        Vector4 pFocus = pLensLocal + (mFocalDistance / dir.z) * dir; 
        Vector3 lensToFilm(
            pFocus.x - pLensLocal.x,
            pFocus.y - pLensLocal.y,
            pFocus.z - pLensLocal.z);
        float lensToFilmDistance2 = lensToFilm.squaredLength();
        float cosTheta = normalize(lensToFilm).z;
        float G = cosTheta * cosTheta / lensToFilmDistance2;
        float lensArea = mLensRadius * mLensRadius * PI; 
        float We = mLensRadius > 0.0f ?
            1.0f / (mFilmArea * lensArea * G) :
            1.0f / (mFilmArea * G);
        return We;
    }

    bool PerspectiveCamera::isDelta() const {
        return (mLensRadius == 0.0f);
    }

    Vector3 PerspectiveCamera::worldToScreen(const Vector3&pWorld,
        const Vector3& pLens) const {
        Vector4 pView = mView * Vector4(pWorld, 1.0f);
        Vector4 pLensLocal = mView * Vector4(pLens, 1.0f);
        Vector3 invalidScreenPoint(0.0f, 0.0f, -1.0f);
        // pWorld is behind camera lens
        if (pView.z < 0.0f) {
            return sInvalidPixel;
        }
        // pLens is not on the lens, invalid input
        if (mLensRadius > 0.0f &&
            pLensLocal.x * pLensLocal.x + pLensLocal.y * pLensLocal.y >
            mLensRadius * mLensRadius) {
            return sInvalidPixel;
        }
        Vector4 dir(pView - pLensLocal);
        // pWorld is on the plane of of lens, it won't be in the screen
        if (fabs(dir.z) < 1e-7f) {
            return sInvalidPixel;
        }
        Vector4 pFocus = pLensLocal + (mFocalDistance / dir.z) * dir;
        Vector4 pNDC = mProj * pFocus;
        pNDC /= pNDC.w;
        float screenX = (pNDC.x + 1.0f) * 0.5f * mFilm->getXResolution();
        float screenY = (1.0f - pNDC.y) * 0.5f * mFilm->getYResolution();
        int xStart, xEnd, yStart, yEnd;
        mFilm->getSampleRange(&xStart, &xEnd, &yStart, &yEnd);
        if (screenX < xStart || screenX > xEnd ||
            screenY < yStart || screenY > yEnd) {
            return sInvalidPixel;
        }
        return Vector3(screenX, screenY, pView.z);
    }

    OrthographicCamera::OrthographicCamera(const Vector3& position, 
        const Quaternion& orientation, float zn, float zf, 
        float filmWidth, Film* film):
        Camera(position, orientation, zn, zf, film),
        mFilmWidth(filmWidth) {
        mFilmHeight = mFilmWidth / mAspectRatio;
        mFilmArea = mFilmWidth * mFilmHeight;
        mPixelArea = mFilmArea *
            mFilm->getInvXResolution() * mFilm->getInvYResolution();
        mFilm->setFilmArea(mFilmArea);
        mProj = matrixOrthoLHD3D(mFilmWidth, mFilmHeight, mZNear, mZFar);
        update();
    }

    float OrthographicCamera::generateRay(const Sample& sample, 
        RayDifferential* ray) const {
        float invXRes = mFilm->getInvXResolution();
        float invYRes = mFilm->getInvYResolution();
        float xNDC = +2.0f * sample.imageX * invXRes - 1.0f;
        float yNDC = -2.0f * sample.imageY * invYRes + 1.0f;
        float dxNDC = +2.0f * (sample.imageX + 1.0f) * invXRes - 1.0f;
        float dyNDC = -2.0f * (sample.imageY + 1.0f) * invYRes + 1.0f;

        // from NDC space to view space
        float xView = 0.5f * mFilmWidth * xNDC;
        float yView = 0.5f * mFilmHeight * yNDC;
        float dxView = 0.5f * mFilmWidth * dxNDC;
        float dyView = 0.5f * mFilmHeight * dyNDC;
        Vector3 pView(xView, yView, 0.0f);
        Vector3 pDxView(dxView, yView, 0.0f);
        Vector3 pDyView(xView, dyView, 0.0f);
        Vector3 viewDir(0.0f, 0.0f, 1.0f);
        // view space to world space
        ray->o = mPosition + mOrientation * pView;
        ray->dxOrigin = mPosition + mOrientation * pDxView;
        ray->dyOrigin = mPosition + mOrientation * pDyView;
        ray->d = ray->dxDir = ray->dyDir = mOrientation * viewDir;
        ray->mint = 0.0f;
        ray->maxt = INFINITY;
        ray->depth = 0;
        ray->hasDifferential = true;
        return 1.0f;
    }

    Vector3 OrthographicCamera::samplePosition(const Sample& sample,
        Vector3* surfaceNormal, float* pdfArea) const {
        float invXRes = mFilm->getInvXResolution();
        float invYRes = mFilm->getInvYResolution();
        float xNDC = +2.0f * sample.imageX * invXRes - 1.0f;
        float yNDC = -2.0f * sample.imageY * invYRes + 1.0f;
        // from NDC space to view space
        float xView = 0.5f * mFilmWidth * xNDC;
        float yView = 0.5f * mFilmHeight * yNDC;
        if (pdfArea) {
            *pdfArea = 1.0f / mFilmArea;
        }
        *surfaceNormal = getLook();
        return mOrientation * Vector3(xView, yView, 0.0f) + mPosition;
    }

    Vector3 OrthographicCamera::sampleDirection(const Sample& sample,
        const Vector3& pCamera, float* We, float* pdfW) const {
        *We = 1.0f / mFilmArea;
        if (pdfW) {
            *pdfW = 1.0f;
        }
        return mOrientation * Vector3::UnitZ;
    }

    float OrthographicCamera::pdfPosition(const Vector3& p) const {
        return 1.0f / mFilmArea;
    }

    float OrthographicCamera::pdfDirection(const Vector3& p,
        const Vector3& wo) const {
        return 0.0f;
    }

    float OrthographicCamera::evalWe(const Vector3& pCamera,
        const Vector3& pWorld) const {
        // OrthographicCamera is not even intersectable, this method
        // should never get called anyway since there won't be particle
        // lands on this type of camera
        return 0.0f;
    }

    bool OrthographicCamera::isDelta() const {
        return true;
    }

    Camera* PerspectiveCameraCreator::create(const ParamSet& params, 
        Film* film) const {
        Vector3 position = params.getVector3("position");
        Quaternion orientation = getQuaternion(params);
        float fov = params.getFloat("fov", 60.0f);
        float zn = params.getFloat("near_plane", 0.1f);
        float zf = params.getFloat("far_plane", 1000.0f);
        float lensRadius = params.getFloat("lens_radius", 0.0f);
        float focalDistance = params.getFloat("focal_distance", 1.0f);
        return new PerspectiveCamera(position, orientation, radians(fov), 
            zn, zf, lensRadius, focalDistance, film);
    }


    Camera* OrthographicCameraCreator::create(const ParamSet& params, 
        Film* film) const {
        Vector3 position = params.getVector3("position");
        Quaternion orientation = getQuaternion(params);
        float zn = params.getFloat("near_plane", 0.1f);
        float zf = params.getFloat("far_plane", 1000.0f);
        float filmWidth = params.getFloat("film_width", 35.0f);
        return new OrthographicCamera(position, orientation, zn, zf, 
            filmWidth, film);
    }
}
