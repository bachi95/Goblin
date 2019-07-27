#include "GoblinSampler.h"
#include "GoblinVector.h"

#include <cassert>

namespace Goblin {
void SampleQuota::clear() {
    n1D.clear();
    n2D.clear();
}
    
size_t SampleQuota::size() const {
    size_t totalSize = 0;
    for (size_t i = 0; i < n1D.size(); ++i) {
        totalSize += n1D[i];
    }
    for (size_t i = 0; i < n2D.size(); ++i) {
        totalSize += 2 * n2D[i];
    }
    return totalSize;
}

SampleIndex SampleQuota::requestOneDQuota(uint32_t samplesNum) {
    int nSample = roundToSquare(samplesNum);
    n1D.push_back(nSample);
    return SampleIndex((uint32_t)n1D.size() - 1, nSample);
}

SampleIndex SampleQuota::requestTwoDQuota(uint32_t samplesNum) {
    int nSample = roundToSquare(samplesNum);
    n2D.push_back(roundToSquare(nSample));
    return SampleIndex((uint32_t)n2D.size() - 1, nSample);
}

void Sample::allocateQuota(const SampleQuota& quota) {
    n1D = quota.n1D;
    n2D = quota.n2D;

    uint32_t patternNum = 0;
    patternNum = (uint32_t)(n1D.size() + n2D.size());
    if (patternNum == 0) {
        u1D = nullptr;
        u2D = nullptr;
        return;
    }
    // allocate a chunk of float* buffer that can host all sample patterns
    u1D = new float*[patternNum];
    u2D = u1D + n1D.size();

    // allocate a chunk of float buffer that can host all 1d/2d values
    float* quotaBuffer = new float[quota.size()];
    for (uint32_t i = 0; i < n1D.size(); ++i) {
        u1D[i] = quotaBuffer;
        quotaBuffer += n1D[i];
    }
    for (uint32_t i = 0; i < n2D.size(); ++i) {
        u2D[i] = quotaBuffer;
        quotaBuffer += 2 * n2D[i];
    }
}


Sampler::Sampler(const SampleRange& sampleRange,
    int samplePerPixel, const SampleQuota& sampleQuota,
    RNG* rng):
    mXStart(sampleRange.xStart), mXEnd(sampleRange.xEnd), 
    mYStart(sampleRange.yStart), mYEnd(sampleRange.yEnd),
    mCurrentX(sampleRange.xStart), mCurrentY(sampleRange.yStart),
    mSampleBuffer(nullptr), mJitter(true),
    mSampleQuota(sampleQuota), mRNG(rng) {
    int root;
    mSamplesPerPixel = roundToSquare(samplePerPixel, &root);
    mXPerPixel = mYPerPixel = root;
}

Sampler::~Sampler() {
    if (mSampleBuffer != nullptr) {
        delete [] mSampleBuffer;
        mSampleBuffer = nullptr;
    }
}

int Sampler::maxSamplesPerRequest() const {
    return mSamplesPerPixel;
}

uint64_t Sampler::maxTotalSamples() const {
    return (uint64_t)mSamplesPerPixel * 
        (uint64_t)(mXEnd - mXStart) * 
        (uint64_t)(mYEnd - mYStart);
}

/* 
    * this need some explanations here...
    * to let all the samples in one pixel cover as much sample range,
    * besides stratified sample values based on number of u1d/u2d
    * per sample pattern (for example, n2D[0] has 25 2d samples
    * they will be placed under 5X5 strata, each strata size 1/5 = 0.2)
    * we further stratified those strata based on samplesPerPixel(
    * continue the previous example, let's say we take 
    * 4 samples per pixel, we stratified the (0.2, 0.2) strata to
    * 2X2 sub strata that each sub strata with size 0.2/2 = 0.1)
    * and fill in the fine chopped float to those sub strata,
    * then for each strata we shuffle its sub strata and assign
    * the value to each samplesPerPixel as the result
    *
    */
int Sampler::requestSamples(Sample* samples) {
    if (mCurrentY == mYEnd) {
        return 0;
    }
    if (mSampleBuffer == nullptr) {
        // 4(imageX, imageY, lensU1, lensU2) + 
        // quota size(extra requested 1/2D samples)
        size_t bufferSize = mSamplesPerPixel * (4 + mSampleQuota.size());
        mSampleBuffer = new float[bufferSize];
    }

    float* imageBuffer = mSampleBuffer;
    float* lensBuffer = mSampleBuffer + 2 * mSamplesPerPixel;
    float* quotaBuffer = mSampleBuffer + 4 * mSamplesPerPixel;
    // 1 pixel is 1 strata, gettting 2d stratified as 
    // mSamplesPerPixel sub strata
    stratifiedUniform2D(imageBuffer, 1);
    stratifiedUniform2D(lensBuffer, 1);
    float* currentOffset = quotaBuffer;
    for (size_t i = 0; i < mSampleQuota.n1D.size(); ++i) {
        stratifiedUniform1D(currentOffset, mSampleQuota.n1D[i]);
        currentOffset += mSampleQuota.n1D[i] * mSamplesPerPixel;
    }
    for (size_t i = 0; i < mSampleQuota.n2D.size(); ++i) {
        stratifiedUniform2D(currentOffset, mSampleQuota.n2D[i]);
        currentOffset += 2 * mSampleQuota.n2D[i] * mSamplesPerPixel;
    }
    // shuffle the above stratified result
    shuffle(lensBuffer, mSamplesPerPixel, 2, mRNG);

    float* shuffleBuffer = quotaBuffer;
    for (size_t i = 0; i < mSampleQuota.n1D.size(); ++i) {
        for (uint32_t j = 0; j < mSampleQuota.n1D[i]; ++j) {
            shuffle(shuffleBuffer, mSamplesPerPixel, 1, mRNG);
            shuffleBuffer += mSamplesPerPixel;
        }
    }
    for (size_t i = 0; i < mSampleQuota.n2D.size(); ++i) {
        for (uint32_t j = 0; j < mSampleQuota.n2D[i]; ++j) {
            shuffle(shuffleBuffer, mSamplesPerPixel, 2, mRNG);
            shuffleBuffer += 2 * mSamplesPerPixel;
        }
    }

    // fill in mSampleBuffer to samples
    for (int i = 0; i < mSamplesPerPixel; ++i) {
        samples[i].imageX = mCurrentX + imageBuffer[2 * i];
        samples[i].imageY = mCurrentY + imageBuffer[2 * i + 1];
    }
    for (int i = 0; i < mSamplesPerPixel; ++i) {
        samples[i].lensU1 = lensBuffer[2 * i];
        samples[i].lensU2 = lensBuffer[2 * i + 1];
    } 
    float* fillinBuffer = quotaBuffer;
    for (size_t i = 0; i < mSampleQuota.n1D.size(); ++i) {
        for (uint32_t j = 0; j < mSampleQuota.n1D[i]; ++j) {
            for (int k = 0; k < mSamplesPerPixel; ++k) {
                samples[k].u1D[i][j] = fillinBuffer[k];
            }
            fillinBuffer += mSamplesPerPixel;
        }
    }
    for (size_t i = 0; i < mSampleQuota.n2D.size(); ++i) {
        for (uint32_t j = 0; j < mSampleQuota.n2D[i]; ++j) {
            for (int k = 0; k < mSamplesPerPixel; ++k) {
                samples[k].u2D[i][2 * j] = fillinBuffer[2 * k];
                samples[k].u2D[i][2 * j + 1] = fillinBuffer[2 * k + 1];
            }
            fillinBuffer += 2 * mSamplesPerPixel;
        }
    }
    // another round of shuffle, so that all the u1D u2D in one sample
    // won't have correlation
    for (int i = 0; i < mSamplesPerPixel; ++i) {
        for (size_t j = 0; j < mSampleQuota.n1D.size(); ++j) {
            shuffle(samples[i].u1D[j], mSampleQuota.n1D[j], 1, mRNG);
        }
        for (size_t j = 0; j < mSampleQuota.n2D.size(); ++j) {
            shuffle(samples[i].u2D[j], mSampleQuota.n2D[j], 2, mRNG);
        }
    }

    //debugOutput(samples);

    if (++mCurrentX == mXEnd) {
        mCurrentX= mXStart;
        mCurrentY++;
    }
    return mSamplesPerPixel;
}

void Sampler::debugOutput(Sample* samples) {
    static bool debugFlip = true;
    if (debugFlip) {
        debugFlip = false;
        // debug output, remove me later
        float* currentOffset = mSampleBuffer;
        std::cout << "image sample \n";
        for (int i = 0; i < mSamplesPerPixel; ++i) {
            float x = currentOffset[2 * i];
            float y = currentOffset[2 * i + 1];
            std::cout << "(" << x << ", " << y << ") ";
        }
        std::cout << "\nlens samples\n";
        currentOffset += 2 * mSamplesPerPixel;
        for (int i = 0; i < mSamplesPerPixel; ++i) {
            float x = currentOffset[2 * i];
            float y = currentOffset[2 * i + 1];
            std::cout << "(" << x << ", " << y << ") ";
        } 

        std::cout << "\nintegrator samples\n";
        currentOffset += 2 * mSamplesPerPixel;

        for (size_t i = 0; i < mSampleQuota.n1D.size(); ++i) {
            for (size_t j = 0; j < mSampleQuota.n1D[i]; ++j) {
                for (int k = 0; k < mSamplesPerPixel; ++k) {
                    std::cout << currentOffset[k] << " ";
                }
                std::cout << std::endl;
                currentOffset += mSamplesPerPixel;
            }
            std::cout << std::endl;
        }

        for (size_t i = 0; i < mSampleQuota.n2D.size(); ++i) {
            for (size_t j = 0; j < mSampleQuota.n2D[i]; ++j) {
                for (int k = 0; k < mSamplesPerPixel; ++k) {
                    int index = 2 * k;
                    std::cout <<"(" << currentOffset[index] << 
                        ", "<< currentOffset[index + 1] << ") ";
                }
                currentOffset += 2 * mSamplesPerPixel;
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        for (int i = 0; i < mSamplesPerPixel; ++i) {
            std::cout << "samples " << i << std::endl;
            std::cout << "n1d\n";
            for (size_t j = 0; j < samples[i].n1D.size(); ++j) {
                for (size_t k = 0; k < samples[i].n1D[j]; ++k) {
                    std::cout << samples[i].u1D[j][k] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << "n2d\n";
            for (size_t j = 0; j < samples[i].n2D.size(); ++j) {
                for (size_t k = 0; k < samples[i].n2D[j]; ++k) {
                    std::cout << "(" << samples[i].u2D[j][2 * k] << 
                        ", " << samples[i].u2D[j][2 * k + 1] << ") ";
                }
                std::cout << std::endl;
            }
        }
    }
}


Sample* Sampler::allocateSampleBuffer(size_t bufferSize) {
    Sample* samples = new Sample[bufferSize];
    for (size_t i = 0; i < bufferSize; ++i) {
        samples[i].allocateQuota(mSampleQuota);
    }
    return samples;
}

void Sampler::stratifiedUniform1D(float* buffer, uint32_t n1D) {
    float strataSize = 1.0f / (float)n1D;
    float subStrataSize = strataSize / mSamplesPerPixel;
    for (uint32_t i = 0; i < n1D; ++i) {
        for (int j = 0; j < mSamplesPerPixel; ++j) {
            float nOffset = mJitter ? j + mRNG->randomFloat() : j + 0.5f;
            int index = i * mSamplesPerPixel + j;
            buffer[index] = i * strataSize + nOffset * subStrataSize;
        }
    }
}

void Sampler::stratifiedUniform2D(float* buffer, uint32_t n2D) {
    int root = (int)sqrtf((float)n2D);
    float strataSize = 1.0f / root;
    float subStrataSize = strataSize / mXPerPixel;
    for (uint32_t n = 0; n < n2D; ++n) {
        int uX = n % root;
        int uY = n / root;
        for (int p = 0; p < mSamplesPerPixel; ++p) {
            int pX = p % mXPerPixel;
            int pY = p / mXPerPixel;
            float xOffset = mJitter ? pX + mRNG->randomFloat() : pX + 0.5f;
            float yOffset = mJitter ? pY + mRNG->randomFloat() : pY + 0.5f;
            int index = 2 * (n * mSamplesPerPixel + pY * mXPerPixel + pX);
            buffer[index] = 
                uX * strataSize + xOffset * subStrataSize; 
            buffer[index + 1 ] = 
                uY * strataSize + yOffset * subStrataSize; 
        }
    }
}

CDF1D::CDF1D(const std::vector<float>& f1D): mFunction(f1D) {
    init();
}

CDF1D::CDF1D(const float* f1D, int n): mFunction(f1D, f1D + n) {
    init();
}

void CDF1D::init() {
    size_t n = mFunction.size();
    mDx = 1.0f / n;
    mCount = (int)n;
    mCDF.resize(n + 1);
    mCDF[0] = 0.0f;
    // accumulate up the integral
    for (size_t i = 1; i < n + 1; ++i) {
        mCDF[i] = mCDF[i - 1] + (mFunction[i - 1] * mDx);
    }
    mIntegral = mCDF[n];
    // normalize CDF with integral
    for (size_t i = 1; i < n + 1; ++i) {
        mCDF[i] /= mIntegral;
    }
}

int CDF1D::sampleDiscrete(float u , float* pdf) {
    std::vector<float>::iterator lowBound;
    lowBound = std::lower_bound(mCDF.begin(), mCDF.end(), u);
    int offset = std::max(0, static_cast<int>(lowBound - mCDF.begin() - 1));
    if (pdf) {
        *pdf = (mFunction[offset] / mIntegral) * mDx;
    }
    return offset;
}

float CDF1D::sampleContinuous(float u, float* pdf, int* index) {
    std::vector<float>::iterator lowBound;
    lowBound = std::lower_bound(mCDF.begin(), mCDF.end(), u);
    int offset = std::max(0, static_cast<int>(lowBound - mCDF.begin() - 1));
    float d = (u - mCDF[offset]) / (mCDF[offset + 1] - mCDF[offset]);
    if (pdf) {
        *pdf = mFunction[offset] / mIntegral;
    }
    if (index) {
        *index = offset;
    }
    return (static_cast<float>(offset) + d) / mFunction.size();
}


CDF2D::CDF2D(const float* f2D, int width, int height) {
    std::vector<float> rowIntegrals;
    for (int i = 0; i < height; ++i) {
        CDF1D* colCDF = new CDF1D(f2D + i * width, width);
        rowIntegrals.push_back(colCDF->mIntegral);
        mConditionalDist.push_back(colCDF);
    }
    mMarginalDist = new CDF1D(rowIntegrals);
}

CDF2D::~CDF2D() {
    for (size_t i = 0; i < mConditionalDist.size(); ++i) {
        delete mConditionalDist[i];
        mConditionalDist[i] = nullptr;
    }
    delete mMarginalDist;
    mMarginalDist = nullptr;
}
    
Vector2 CDF2D::sampleContinuous(float u1, float u2, float* pdf) {
    // first pick up the row based on marginal pdf alone rows
    float pdfRow;
    int row;
    float v = mMarginalDist->sampleContinuous(u2, &pdfRow, &row);
    // the conditional pdf under the condition that we pick row from above
    float pdfCol;
    float u = mConditionalDist[row]->sampleContinuous(u1, &pdfCol);
    if (pdf) {
        *pdf = pdfRow * pdfCol;
    }
    return Vector2(u, v);
}

float CDF2D::pdf(float u, float v) {
    int row = clamp(floorInt(mMarginalDist->mCount * v), 
        0, mMarginalDist->mCount - 1);
    int col = clamp(floorInt(mConditionalDist[row]->mCount * u), 
        0, mConditionalDist[row]->mCount - 1);
    float integral = 
        mMarginalDist->mIntegral * mConditionalDist[row]->mIntegral;
    if (integral == 0.0f) {
        return 0.0f;
    }
    float pdf = mMarginalDist->mFunction[row] * 
        mConditionalDist[row]->mFunction[col] / integral;
    return pdf;
}

/*
    * for a isosceles right triangle with area 1/2 
    * (derivation is the same for other case)
    * pdf(u, v) = 1/(1/2) = 2 (pdf with respect to area)
    * marginal pdf(u) = integrate pdf(u, v) over 0 to 1 -u (0 <= v <= 1 - u)
    * = 2(1 - u) 
    * conditional pdf(u|v) = pdf(u, v)/pdf(u) = 1 / (1 - u)
    * cdf(u) = integratge pdf(u) over 0 to u = 2u - u^2
    * cdf(v) = integrate pdf(v|u) over 0 to v = v / (1 - u)
    * inverse method:
    * u1 = 2u - u^2 -> u = 1 - sqrt(1-u1) can be simplified as 1 -sqrt(u1)
    * u2 = v / (1-u) -> v = sqrt(u1) * u2
    */
void uniformSampleTriangle(float u1, float u2, float* u, float* v) {
    float u1root = sqrtf(u1);
    *u = 1.0f - u1root;
    *v = u1root * u2;
}
    
/*
    * for cone 0 <= theta <= thetaMax, 0 <= phi <= 2pi
    * the solid angle extended by the cone is:
    * integrate dw = sin(theta)d(theta)d(phi) over theta and phi
    * = 2pi(1 - cos(thetaMax))
    * pdf(w) = 1 / 2pi(1 - cos(thetaMax))
    * pdf(theta, phi) = sin(theta) / 2pi(1 - cos(thetaMx))
    * marginal pdf(theta) = integrate pdf(theta, phi) over 0 to 2pi
    * = sin(theta) / (1 - cos(thetaMax))
    * conditional pdf(phi|theta) = pdf(theta, phi) / pdf(theta) = 1 / 2pi
    * cdf(theta) = integrate sin(theta) / (1 - cos(thetaMax)) over 0 to theta
    * cdf(theta) = (1 - cos(theta)) / (1 - cos(thetaMax))
    * cdf(phi) = integrate 1 / 2pi over 0 to phi = phi / 2pi
    * inverse method:
    * u1 = (1 - cos(theta)) / (1 - cos(thetaMax)) ->
    * cos(theta) = 1 - u1 + u1 * cos(thetaMax)
    * u2 = phi / 2pi -> phi = 2pi * u2
    * transofrm it back to spherical coordinate:
    * sin(theta) = sqrt(1 - cos^2(theta))
    * x = cos(phi) * sin(theta)
    * y = sin(phi) * sin(theta)
    * z = cos(theta)
    */
Vector3 uniformSampleCone(float u1, float u2, float cosThetaMax) {
    float cosTheta = 1.0f - u1 + u1 * cosThetaMax;
    float sinTheta = sqrtf(std::max(0.0f, 1.0f - cosTheta * cosTheta));
    float phi = TWO_PI * u2;
    float z = cosTheta;
    float x = sinTheta * cos(phi);
    float y = sinTheta * sin(phi);
    return Vector3(x, y, z);
}

Vector3 uniformSampleCone(float u1, float u2, float cosThetaMax,
    const Vector3& x, const Vector3& y, const Vector3& z) {
    float cosTheta = 1.0f - u1 + u1 * cosThetaMax;
    float sinTheta = sqrtf(std::max(0.0f, 1.0f - cosTheta * cosTheta));
    float phi = TWO_PI * u2;
    return x * sinTheta * cos(phi) + 
        y * sinTheta * sin(phi) +
        z * cosTheta;
}

/*
    * uniform means integrate pdf(w) over sphere = 1
    * pdf(w) = 1 / 4pi, since dw = sin(theta)d(theta)d(phi)
    * pdf(theta, phi) = sin(theta) / 4pi
    * marginal pdf(theta) = integrate pdf(theata, phi) over 0 to 2pi
    * = sin(theta) / 2
    * conditional pdf(phi|theta) = pdf(theta, phi) / pdf(theta) = 1 / 2pi
    * cdf(theta) = integrate sin(theta) / 2 over 0 to theta ->
    * cdf(theta) = (1 - cos(theta)) / 2
    * cdf(phi) = integrate 1 / 2pi over 0 to phi = phi / 2pi
    * inverse method:
    * u1 = (1 - cos(theta)) / 2  -> theta = arcos(1 - 2 * u1) 
    * u2 = phi / 2pi -> phi = 2pi * u2
    * transform it back to spherical coordinate:
    * cos(theta) = 1 - 2 * u1
    * sin(theta) = sqrt(1 - cos^2(theta))
    * x = cos(phi) * sin(theta)
    * y = sin(phi) * sin(theta)
    * z = cos(theta)
    */
Vector3 uniformSampleSphere(float u1, float u2) {
    float z = 1.0f - 2.0f * u1;
    float sinTheta = sqrtf(std::max(0.0f, 1.0f - z * z));
    float phi = TWO_PI * u2;
    float x = sinTheta * cos(phi);
    float y = sinTheta * sin(phi);
    return Vector3(x, y, z);
}

/*
    * uniform means integrate pdf(w) over hemisphere = 1
    * pdf(w) = 1 / 2pi, since dw = sin(theta)d(theta)d(phi)
    * pdf(theta, phi) = sin(theta) / 2pi
    * marginal pdf(theta) = integrate pdf(theata, phi) over 0 to 2pi
    * = sin(theta)
    * conditional pdf(phi|theta) = pdf(theta, phi) / pdf(theta) = 1 / 2pi
    * cdf(theta) = integrate sin(theta) over 0 to theta = 1 - cos(theta)
    * cdf(phi) = integrate 1 / 2pi over 0 to phi = phi / 2pi
    * inverse method:
    * u1 = 1 - cos(theta) -> theta = arcos(1 - u1) can be simplified as
    * theta = arcos(u1) (since u1 is uniform distribution)
    * phi = 2pi * u2
    * transform it back to spherical coordinate:
    * cos(theta) = u1
    * x = cos(phi) * sin(theta) = cos(2pi * u2) * sqrt(1 - u1^2)
    * y = sin(phi) * sin(theta) = sin(2pi * u2) * sqrt(1 - u1^2)
    * z = cos(theta) = u1
    */
Vector3 uniformSampleHemisphere(float u1, float u2) {
    float z = u1;
    float sinTheta = sqrtf(std::max(0.0f, 1.0f - u1 * u1));
    float phi = TWO_PI * u2;
    float x = sinTheta * cos(phi);
    float y = sinTheta * sin(phi);
    return Vector3(x, y, z);
}

/*
    * since the pdf is proportional to cos(theta)
    * let pdf(w) = c * cos(theta)
    * pdf(theta, phi) = c * sin(theta) * cos(tehta) = (c/2)  * sin(2 * theta)
    * integrate c/2 * sin(2*theta) over theta and phi = 1 ->
    * c = 1 / pi -> pdf(w) = cos(theta) / pi ->
    * pdf(theta, phi) =  sin(2 * theta) / 2pi
    * marginal pdf(theta) = integrage pdf(theta, phi) over 0 to 2pi
    * = sin(2 * theta)
    * conditional pdf(phi|theta) = 1 / 2pi
    * cdf(theta) = integrate sin(theta) * cos(theta) over 0 to theta ->
    * cdf(theta) = (1 - cos(2 * theta)) / 2 = sin^2(theta)
    * cdf(phi) = integrate 1 / 2pi over 0 to phi = phi / 2pi
    * inverse method:
    * u1 = sin^2(theta) -> theta = arcsin(sqrt(u1))
    * u2 = phi / 2pi -> phi = 2pi * u2
    * transform it back to spherical coordinate:
    * sinTheta = sqrt(u1)
    * cosTheta = sqrt(1 - sinTheta * sinTheta) = sqrt(1 - u1)
    * x = cos(phi) * sin(theta)
    * y = sin(phi) * sin(theta)
    * z = cos(theta)
    */
Vector3 cosineSampleHemisphere(float u1, float u2) {
    float sinTheta = sqrtf(u1);
    float cosTheta = sqrtf(std::max(0.0f, 1.0f - u1));
    float phi = TWO_PI * u2;
    float z = cosTheta;
    float x = sinTheta * cos(phi);
    float y = sinTheta * sin(phi);
    return Vector3(x, y, z);
}

/* 
    * based on Shirley. P, Chiu. K 1997
    * "A low distortion map between disk and square"
    * mapping x y coordinate from square to r, theta in disk
    * r = x, theta = y / x
    */
Vector2 uniformSampleDisk(float u1, float u2) {
    float r, theta;
    // transform [0, 1] sample to [-1, 1]
    float x = 2.0f * u1 - 1.0f;
    float y = 2.0f * u2 - 1.0f;

    if (x + y > 0) {
        if (x > y) {
            // right quarter of square
            r = x;
            theta = 0.25f * PI * (y / x);
        } else {
            // up quarter of square
            r = y;
            // (-pi / 4) * (x / y) + pi / 2
            theta = 0.25f * PI * (2.0f - x / y);
        }
    } else {
        if (x < y) {
            // left quarter of square
            r = -x;
            // (pi / 4) * (y / x) + pi
            theta = 0.25f * PI * (4.0f + y / x);
        } else {
            // down quarter of square
            r = -y;
            // x and y may both be 0
            if (y != 0.0f) {
                // (-pi / 4) * (x / y) + (3pi / 2)
                theta = 0.25f * PI * (6.0f - x / y);
            } else {
                theta = 0.0f;
            }
        }
    }

    return Vector2(r * cos(theta), r * sin(theta));
}

/*
    * let falloff be a
    * since the pdf is proportional to exp(-a * (x * x + y * y))
    * let pdf(x, y) = c * exp(-a * (x * x + y * y))
    * convert to spherical coordinate:
    * pdf(r, theta) = c * r * exp(-a * r * r)
    * integrate c * r * exp(-a * r * r) over r (from 0 to inf) and theta = 1 ->
    * c = a / pi -> pdf(r, theta) = a / pi * r * exp(-a * r * r) ->
    * pdf(x, y) = a / pi * exp(-a * (x * x + y * y))
    * marginal pdf(r) = integrage pdf((r, theta) over 0 to 2pi
    * = 2 * a * r * exp(-a * r * r)
    * conditional pdf(theta|r) = 1 / 2pi
    * cdf(r) = integrate 2 * a * r * exp(-a * r * r) over 0 to r ->
    * cdf(r) = 1 - exp(-a * r * r)
    * cdf(theta) = integrate 1 / 2pi over 0 to theta = theta / 2pi
    * inverse method:
    * u1 = 1 - exp(-a * r * r) -> r = sqrt(ln(1 - u1) / -a) ->
    * r = sqrt(ln(u1) / -a) since u1 is 0-1 uniform distribution
    * u2 = theta / 2pi -> theta = 2pi * u2
    * transform it back from spherical coordinate:
    * x = r * cos(theta)
    * y = r * sin(theta)
    */
Vector2 gaussianSample2D(float u1, float u2, float falloff) {
    float r = sqrtf(log(u1) / -falloff);
    float theta = TWO_PI * u2;
    return Vector2(r * cos(theta), r * sin(theta));
}

/*
    * similar to the above gaussianSample2D, but the integrate range of r
    * from 0 to Rmax, which makes us able to sample a disc with radius Rmax in
    * gaussian distribution
    */
Vector2 gaussianSample2D(float u1, float u2, float falloff, float Rmax) {
    float r = sqrtf(log(1.0f - u1 * (1.0f - exp(-falloff * Rmax * Rmax))) /
        -falloff);
    float theta = TWO_PI * u2;
    return Vector2(r * cos(theta), r * sin(theta));
}

float gaussianSample2DPdf(const Vector3& pCenter, 
    const Vector3& pSample, const Vector3& N, float falloff) {
    Vector3 d = pSample - pCenter;
    Vector3 projected = d - N * dot(d, N);
    return INV_PI * falloff * exp(-falloff * squaredLength(projected));
}


float gaussianSample2DPdf(const Vector3& pCenter, 
    const Vector3& pSample, const Vector3& N, float falloff, float Rmax) {
    return gaussianSample2DPdf(pCenter, pSample, N, falloff) /
        (1.0f - exp(-falloff * Rmax * Rmax));
}
    
PermutedHalton::PermutedHalton(size_t dimension, RNG* rng) {
    getPrimes(dimension, mPrimes);
    mTableIndexes.resize(dimension);
    size_t tableSize = 0;
    for (size_t i = 0; i < mPrimes.size(); ++i) {
        tableSize += mPrimes[i];
    }
    mPermutedTable.resize(tableSize);
    size_t offset = 0;
    for (size_t i = 0; i < mPrimes.size(); ++i) {
        mTableIndexes[i] = offset;
        uint32_t p = mPrimes[i];
        for (uint32_t j = 0; j < p; ++j) {
            mPermutedTable[offset + j] = j;
        }
        shuffle(&mPermutedTable[offset], p, 1, rng);
        offset += p;
    }
}

void PermutedHalton::sample(Sample* s, int pixelX, int pixelY,
    uint64_t n, RNG*rng) const {
    size_t dimension = mPrimes.size();
    s->imageX = dimension < 1 ? pixelX + rng->randomFloat() :
        pixelX + permutedRadicalInverse(n, mPrimes[0],
        &mPermutedTable[mTableIndexes[0]]);
    s->imageY = dimension < 2 ? pixelY + rng->randomFloat() :
        pixelY + permutedRadicalInverse(n, mPrimes[1],
        &mPermutedTable[mTableIndexes[1]]);
    s->lensU1 = dimension < 3 ? rng->randomFloat() :
        permutedRadicalInverse(n, mPrimes[2],
        &mPermutedTable[mTableIndexes[2]]);
    s->lensU2 = dimension < 4 ? rng->randomFloat() :
        permutedRadicalInverse(n, mPrimes[3],
        &mPermutedTable[mTableIndexes[3]]);
    size_t currentDimIndex = 4;
    for (size_t i = 0; i < s->n1D.size(); ++i) {
        for (size_t j = 0; j < s->n1D[i]; ++j) {
            s->u1D[i][j] = currentDimIndex >= dimension ?
                rng->randomFloat() :
                permutedRadicalInverse(n, mPrimes[currentDimIndex],
                &mPermutedTable[mTableIndexes[currentDimIndex]]);
            currentDimIndex++;
        }
    }
    for (size_t i = 0; i <s->n2D.size(); ++i) {
        for (size_t j = 0; j < s->n2D[i]; ++j) {
            s->u2D[i][2 *j] = currentDimIndex >= dimension ?
                rng->randomFloat() :
                permutedRadicalInverse(n, mPrimes[currentDimIndex],
                &mPermutedTable[mTableIndexes[currentDimIndex]]);
            currentDimIndex++;
            s->u2D[i][2 *j + 1] = currentDimIndex >= dimension ?
                rng->randomFloat() :
                permutedRadicalInverse(n, mPrimes[currentDimIndex],
                &mPermutedTable[mTableIndexes[currentDimIndex]]);
            currentDimIndex++;
        }
    }
}

void PermutedHalton::sample(Sample* s, uint64_t n, RNG*rng) const {
    s->imageX = s->imageY = s->lensU1 = s->lensU2 = 0.0f;
    size_t dimension = mPrimes.size();
    size_t currentDimIndex = 0;
    for (size_t i = 0; i < s->n1D.size(); ++i) {
        for (size_t j = 0; j < s->n1D[i]; ++j) {
            s->u1D[i][j] = currentDimIndex >= dimension ?
                rng->randomFloat() :
                permutedRadicalInverse(n, mPrimes[currentDimIndex],
                &mPermutedTable[mTableIndexes[currentDimIndex]]);
            currentDimIndex++;
        }
    }
    for (size_t i = 0; i <s->n2D.size(); ++i) {
        for (size_t j = 0; j < s->n2D[i]; ++j) {
            s->u2D[i][2 *j] = currentDimIndex >= dimension ?
                rng->randomFloat() :
                permutedRadicalInverse(n, mPrimes[currentDimIndex],
                &mPermutedTable[mTableIndexes[currentDimIndex]]);
            currentDimIndex++;
            s->u2D[i][2 *j + 1] = currentDimIndex >= dimension ?
                rng->randomFloat() :
                permutedRadicalInverse(n, mPrimes[currentDimIndex],
                &mPermutedTable[mTableIndexes[currentDimIndex]]);
            currentDimIndex++;
        }
    }
}
}
