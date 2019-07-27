#ifndef GOBLIN_SAMPLER_H
#define GOBLIN_SAMPLER_H

#include "GoblinUtils.h"

namespace Goblin {
class Vector2;
class Vector3;
class CDF2D;

// use for book keeping where and how many u1d/u2d can be
// retrieved from Sample
struct SampleIndex {
    SampleIndex() {};
    SampleIndex(uint32_t o, uint32_t n):
        offset(o), sampleNum(n) {}
    uint32_t offset;
    uint32_t sampleNum;
};

class SampleQuota {
public:
    void clear();
    SampleIndex requestOneDQuota(uint32_t samplesNum);
    SampleIndex requestTwoDQuota(uint32_t samplesNum);
    size_t size() const;
    // how many "float" needs for Sample based on this SampleQuota
    // (includes the 4 extra dimensions for pixel(2) and lens(2)
    size_t getDimension() const { return size() + 4; }

    std::vector<uint32_t> n1D;
    std::vector<uint32_t> n2D;
};


class Sample {
public:
    Sample();
    ~Sample();
    void allocateQuota(const SampleQuota& quota);
    // store film sample in image space (not NDC space)
    float imageX, imageY;
    // used to sample lens for DOF
    float lensU1, lensU2;
    std::vector<uint32_t> n1D;
    std::vector<uint32_t> n2D;
    float** u1D;
    float** u2D;
};

inline Sample::Sample():
    imageX(0.0f), imageY(0.0f), lensU1(0.0f), lensU2(0.0f), 
    u1D(nullptr), u2D(nullptr) {
    n1D.clear();
    n2D.clear();
}

inline Sample::~Sample() {
    // release u1D/u2D
    if (u1D != nullptr) {
        // release the float buffer that u1D/u2D point to
        if (u1D[0] != nullptr) {
            delete [] u1D[0];
            u1D[0] = nullptr;
        }
        delete [] u1D;
        u1D = nullptr;
    }
}

struct SampleRange {
    SampleRange(): xStart(0), xEnd(0), yStart(0), yEnd(0) {}

    SampleRange(int xs, int xe, int ys, int ye):
        xStart(xs), xEnd(xe), yStart(ys), yEnd(ye) {}

    int xStart;
    int xEnd;
    int yStart;
    int yEnd;
};

class Sampler {
public:
    Sampler(const SampleRange& sampleRange, 
        int samplePerPixel, const SampleQuota& sampleQuota,
        RNG* rng);
    ~Sampler();
    int maxSamplesPerRequest() const;
    uint64_t maxTotalSamples() const;
    int requestSamples(Sample* samples);

    Sample* allocateSampleBuffer(size_t bufferSize);
private:
    void stratifiedUniform1D(float* buffer, uint32_t n1D);
    void stratifiedUniform2D(float* buffer, uint32_t n2D);

    void debugOutput(Sample* samples);
private:
    int mXStart, mXEnd;
    int mYStart, mYEnd;
    int mCurrentX, mCurrentY;
    int mXPerPixel, mYPerPixel;
    int mSamplesPerPixel;
    float* mSampleBuffer;
    bool mJitter;
    SampleQuota mSampleQuota;
    RNG* mRNG;
};

/*
    * Cumulative Distribution Function 1D
    * feed in a 1d function in vector form
    * construct its CDF range from 0 to 1
    * and offer interface for sampling based on 
    * the CDF( ie. higher function value, higher pdf )
    * this will be used for sampling light sources
    * based on their power distribution and smapling 
    * geometry from geometries based on their area distribution
    * f(1)=1 f(2)=3 -> dx=1/2 integral=(1+3)*dx=2
    * CDF(0)=0 CDF(1)=1*dx/integral=0.25 CDF(2)=(1+3)*dx/integral=1
    */
class CDF1D {
public:
    CDF1D(const std::vector<float>& f1D);
    CDF1D(const float* f1D, int n);
    int sampleDiscrete(float u, float* pdf = nullptr);
    float sampleContinuous(float u, float* pdf = nullptr, int *index = nullptr);
private:
    void init();
private:
    std::vector<float> mFunction;
    std::vector<float> mCDF;
    float mIntegral;
    float mDx;
    int mCount;
    friend class CDF2D;
};

class CDF2D {
public:
    CDF2D(const float* f2D, int width, int height);
    ~CDF2D();
    Vector2 sampleContinuous(float u1, float u2, float* pdf = nullptr);
    float pdf(float u, float v);
private:
    CDF1D* mMarginalDist;
    std::vector<CDF1D*> mConditionalDist;
};

template<typename T>
void shuffle(T* buffer, uint32_t num, uint32_t dim, RNG* rng) {
    for (uint32_t n = 0; n < num; ++n) {
        size_t toShuffle = rng->randomUInt() % num;
        for (uint32_t d = 0; d < dim; ++d) {
			std::swap(buffer[n * dim + d], buffer[toShuffle * dim + d]);
        }
    }
}

/*
    * u1, u2 are 2 [0, 1) sample value in uniform distribution
    * u, v are barycentric coordinates for two vertices in triangle 
    */
void uniformSampleTriangle(float u1, float u2, float* u, float* v);

Vector3 uniformSampleCone(float u1, float u2, float cosThetaMax);

Vector3 uniformSampleCone(float u1, float u2, float cosThetaMax,
    const Vector3& x, const Vector3& y, const Vector3& z);

inline float uniformConePdf(float cosThetaMax) {
    return 1.0f / (TWO_PI * (1.0f - cosThetaMax));
}

Vector3 uniformSampleSphere(float u1, float u2);

inline float uniformSpherePdf() {
    // 1 / 4pi
    return  0.5f * INV_TWOPI;
}
Vector3 uniformSampleHemisphere(float u1, float u2);

inline float uniformHemispherePdf() {
    return INV_TWOPI;
}

Vector3 cosineSampleHemisphere(float u1, float u2);

inline float cosineHemispherePdf(float cosTheta) {
    return  cosTheta * INV_PI;
}

Vector2 uniformSampleDisk(float u1, float u2);

Vector2 gaussianSample2D(float u1, float u2, float falloff);

inline float gaussianSample2DPdf(float x, float y, float falloff) {
    return INV_PI * falloff * exp(-falloff * (x * x + y * y));
}

Vector2 gaussianSample2D(float u1, float u2, float falloff, float Rmax);

inline float gaussianSample2DPdf(float x, float y, float falloff, 
    float Rmax) {
    return gaussianSample2DPdf(x, y, falloff) / 
        (1.0f - exp(-falloff * Rmax * Rmax));
}

// project pSample on the plane form by pCenter and N, then calculate
// the corresponding gaussian 2D pdf 
float gaussianSample2DPdf(const Vector3& pCenter, 
    const Vector3& pSample, const Vector3& N, float falloff);

float gaussianSample2DPdf(const Vector3& pCenter, 
    const Vector3& pSample, const Vector3& N, float falloff, float Rmax);

/*
    * integrate c * exp(-falloff * x) from 0 to inf = 1 ->
    * c = falloff -> pdf(x) = falloff * exp(-falloff * x)
    * cdf(x) = 1 - exp(-falloff * x) -> u = 1 - exp(-falloff * x) ->
    * x = -ln(1 - u) / falloff -> simplified to x = -ln(u) / falloff
    * since u are [0, 1) uniform distribution
    */
inline float exponentialSample(float u, float falloff) {
    return -log(u) / falloff;
}

inline float exponentialPdf(float x, float falloff) {
    return falloff * exp(-falloff * x);
}

/*
    * similar to above exponentialSample but instead of drawing a sample
    * range between 0 to infinite, this one draw a sample t between a and b
    * that its pdf is propotional to exp(-sigma * (t - a))
    * integrate c * exp(-sigma *(t - 1) from a to b = 1->
    * c * exp(sigma * a) * (exp(-sigma * b) - exp(-sigma * a) / (-sigma) = 1
    * c = -sigma / (exp(sigma * (a - b)) - 1)
    * pdf(t) = -sigma * exp(sigma * a) * exp(-sigma * t) /
    *          (exp(sigma * (a - b) - 1) =
    *          sigma / (exp(sigma * (t - a)) - exp(sigma * (t - b)))
    * cdf(t) = integrate pdf from a to t =
    *          (exp(sigma * (a - t)) - 1) / (exp(sigma * (a - b)) - 1)
    * inverse method:
    * u = (exp(sigma * (a - t)) - 1) / (exp(sigma * (a - b)) - 1)
    * u * (exp(sigma * (a - b)) - 1) + 1 = exp(sigma * (a - t))
    * sigma * (a - t) = log(u * (exp(sigma * (a - b)) - 1) + 1)
    * t = a - log(1 - u * (1 - exp(sigma * (a - b)))) / sigma
    */
inline float exponentialSample(float u, float sigma, float a, float b) {
    return a - log(1.0f - u * (1.0f - exp(sigma * (a - b)))) / sigma;
}

inline float exponentialPdf(float t, float sigma, float a , float b) {
    return sigma / (exp(sigma * (t - a)) - exp(sigma * (t - b)));
}

/*
    * see Kulla. C, Farjardo. M 2012
    * "Importance Sampling Techniques for Path Tracing in Participating Media"
    * for detail reference
    * The basic idea is we try to sample a 1d t value alone the ray that
    * the point distribution is proportional to its inverse square distance
    * to light:
    *                        *
    *                       /|\
    *                      / |  \
    *                     /  |    \
    *                    /   |D     \
    *                   /    |        \
    * -----------------|-----|---|-----|---------
    *                  a         t     b
    * pdf(t) = c / (D^2 + t^2)
    * cdf(t) = integrate pdf from a to b =
    *          (c / D) * (arctan(a / D) - arctan(b / D)) =
    *          (c / D) * (thetaA - thetaB) = 1 ->
    * c = D / (thetaA - thetaB) ->
    * pdf(t) = D / ((thetaB - thetaA) * (D^2 + t^2))
    * inverse method:
    * cdf(t) = (1 / (thetaB - thetaA)) * (arctan(t / D) - thetaA) = u
    * D * tan((thetaB - thetaA) * u + thetaA) = t ->
    * t = D *tan((1 - u) * thetaA + u * thetaB)
    */
inline float equiAngularSample(float u, float D,
    float thetaA, float thetaB) {
    return D * tan((1 - u) * thetaA + u * thetaB);
}

inline float equiAngularPdf(float t, float D, float thetaA, float thetaB) {
    return D / ((thetaB - thetaA) * (D * D + t * t));
}

// determine weight for multi-importance sampling
inline float powerHeuristic(float nA, float pdfA, float nB, float pdfB) {
    float A = nA * pdfA;
    float B = nB * pdfB;
    return A * A / (A * A + B * B);
}

inline float radicalInverse(uint64_t N, uint32_t base) {
    float invBase = 1.0f / base;
    float invBi = invBase;
    float result = 0.0f;
    while (N > 0) {
        uint64_t di = N % base;
        result += di * invBi;
        N /= base;
        invBi *= invBase;
    }
    return clamp(result, 0.0f, 1.0f);
}

inline float permutedRadicalInverse(uint64_t N, uint32_t base,
    const uint32_t* const permutedTable) {
    float invBase = 1.0f / base;
    float invBi = invBase;
    float result = 0.0f;
    while (N > 0) {
        uint64_t di = permutedTable[N % base];
        result += di * invBi;
        N /= base;
        invBi *= invBase;
    }
    // if 0 got permuted to somewhere else than permutedTable[0],
    // we need to take the infinite series of 0 trailing after last di
    // into account. the infinite series will be:
    // permutedTable[0] * invBi * (1 + 1/base + 1/base^2 + 1/base^3...) =
    // permutedTable[0] * invBi * (1 / (1 - 1/base)) =
    // permutedTable[0] * invBi * base / (base - 1)
    result += permutedTable[0] * invBi * base / (base - 1.0f);
    return clamp(result, 0.0f, 1.0f);
}

class PermutedHalton {
public:
    PermutedHalton(size_t dimension, RNG* rng);

    // fill out a Sample based on the sequence id n
    // if the sample dimension is higher than this PermutedHalton,
    // use the feed in RNG to fill up rest of the dimension
    void sample(Sample* s, int pixelX, int pixelY,
        uint64_t n, RNG* rng) const;

    // same as above method but not fill in image/lens pixel related sample
    // can be used for sample that is not emitted from camera
    void sample(Sample* s, uint64_t n, RNG* rng) const;

private:
    std::vector<uint32_t> mPermutedTable;
	std::vector<uint32_t> mPrimes;
    std::vector<size_t> mTableIndexes;
};
}

#endif //GOBLIN_SAMPLER_H


