#ifndef GOBLIN_SAMPLER_H
#define GOBLIN_SAMPLER_H

#include "GoblinUtils.h"

namespace Goblin {
    class Vector2;
    class Vector3;

    struct SampleQuota {
        void clear();
        size_t size() const;
        vector<uint32_t> n1D;
        vector<uint32_t> n2D;
    };

    // use for book keeping where and how many u1d/u2d can be
    // retrieved from Sample
    struct SampleIndex {
        SampleIndex() {};
        SampleIndex(uint32_t o, uint32_t n):
            offset(o), sampleNum(n) {}
        uint32_t offset;
        uint32_t sampleNum;
    };

    class Sample {
    public:
        Sample();
        ~Sample();
        void allocateQuota(const SampleQuota& quota);
        // store film sample in image space (not NDC space)
        float imageX, imageY;

        vector<uint32_t> n1D;
        vector<uint32_t> n2D;
        float** u1D;
        float** u2D;
    };

    inline Sample::Sample():
        imageX(0), imageY(0), u1D(NULL), u2D(NULL) {
        n1D.clear();
        n2D.clear();
    }

    inline Sample::~Sample() {
        // release u1D/u2D
        if(u1D != NULL) {
            // release the float buffer that u1D/u2D point to
            if(u1D[0] != NULL) {
                delete [] u1D[0];
                u1D[0] = NULL;
            }
            delete [] u1D;
            u1D = NULL;
        }
    }

    class Sampler {
    public:
        Sampler(int xStart, int xEnd, int yStart, int yEnd, 
            int samplePerPixel);
        ~Sampler();
        int maxSamplesPerRequest() const;
        int maxTotalSamples() const;
        int requestSamples(Sample* samples);

        SampleIndex requestOneDQuota(uint32_t samplesNum);
        SampleIndex requestTwoDQuota(uint32_t samplesNum);
        Sample* allocateSampleBuffer(size_t bufferSize);
    private:
        void stratifiedUniform1D(float* buffer, uint32_t n1D);
        void stratifiedUniform2D(float* buffer, uint32_t n2D);
        void shuffle(float* buffer, uint32_t num, uint32_t dim);

        void debugOutput(Sample* samples);
    private:
        int mXStart, mYStart;
        int mXEnd, mYEnd;
        int mCurrentX, mCurrentY;
        int mXPerPixel, mYPerPixel;
        int mSamplesPerPixel;
        float* mSampleBuffer;
        bool mJitter;
        SampleQuota mSampleQuota;
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
        CDF1D(const vector<float>& f1D);
        int sampleDiscrete(float u, float* pdf = NULL);
    private:
        vector<float> mFunction;
        vector<float> mCDF;
        float mIntegral;
        float mDx;
    };


    /*
     * u1, u2 are 2 [0, 1) sample value in uniform distribution
     * u, v are barycentric coordinates for two vertices in triangle 
     */
    void uniformSampleTriangle(float u1, float u2, float* u, float* v);

    Vector3 uniformSampleCone(float u1, float u2, float cosThetaMax);

    Vector3 uniformSampleSphere(float u1, float u2);

    Vector3 uniformSampleHemisphere(float u1, float u2);

    Vector3 cosineSampleHemisphere(float u1, float u2);

    Vector2 uniformSampleDisk(float u1, float u2);

    inline float uniformConePdf(float cosThetaMax) {
        return 1.0f / (TWO_PI * (1.0f - cosThetaMax));
    }

    inline float uniformSpherePdf() {
        // 1 / 4pi
        return  0.5f * INV_TWOPI;
    }

    inline float uniformHemispherePdf() {
        return INV_TWOPI;
    }

    inline float cosineHemispherePdf(float cosTheta) {
        return  cosTheta * INV_PI;
    }

    // determine weight for multi-importance sampling
    inline float powerHeuristic(float nA, float pdfA, float nB, float pdfB) {
        float A = nA * pdfA;
        float B = nB * pdfB;
        return A * A / (A * A + B * B);
    }
}

#endif //GOBLIN_SAMPLER_H


