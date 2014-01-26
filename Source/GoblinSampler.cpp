#include "GoblinSampler.h"

namespace Goblin {
    void SampleQuota::clear() {
        n1D.clear();
        n2D.clear();
    }
    
    size_t SampleQuota::size() const {
        size_t totalSize = 0;
        for(size_t i = 0; i < n1D.size(); ++i) {
            totalSize += n1D[i];
        }
        for(size_t i = 0; i < n2D.size(); ++i) {
            totalSize += 2 * n2D[i];
        }
        return totalSize;
    }

    void Sample::allocateQuota(const SampleQuota& quota) {
        n1D = quota.n1D;
        n2D = quota.n2D;

        uint32_t patternNum = 0;
        patternNum = n1D.size() + n2D.size();
        if(patternNum == 0) {
            u1D = NULL;
            u2D = NULL;
            return;
        }
        // allocate a chunk of float* buffer that can host all sample patterns
        u1D = new float*[patternNum];
        u2D = u1D + n1D.size();

        // allocate a chunk of float buffer that can host all 1d/2d values
        float* quotaBuffer = new float[quota.size()];
        for(uint32_t i = 0; i < n1D.size(); ++i) {
            u1D[i] = quotaBuffer;
            quotaBuffer += n1D[i];
        }
        for(uint32_t i = 0; i < n2D.size(); ++i) {
            u2D[i] = quotaBuffer;
            quotaBuffer += 2 * n2D[i];
        }
    }


    Sampler::Sampler(int xStart, int xEnd, int yStart, int yEnd,
        int samplePerPixel):
        mXStart(xStart), mXEnd(xEnd), 
        mYStart(yStart), mYEnd(yEnd),
        mCurrentX(xStart), mCurrentY(yStart),
        mSampleBuffer(NULL), mJitter(false) {

        int root;
        mSamplesPerPixel = roundToSquare(samplePerPixel, &root);
        mXPerPixel = mYPerPixel = root;
        std::cout << "round sample to " << mSamplesPerPixel << 
            " per pixels\n";
        std::cout << "(x, y) sample = (" << 
            mXPerPixel << ", " << mYPerPixel << ")\n";
        mSampleQuota.clear();
    }

    Sampler::~Sampler() {
        if(mSampleBuffer != NULL) {
            delete [] mSampleBuffer;
            mSampleBuffer = NULL;
        }
    }

    int Sampler::maxSamplesPerRequest() const {
        return mSamplesPerPixel;
    }

    int Sampler::maxTotalSamples() const {
        return mSamplesPerPixel * (mXEnd - mXStart) * (mYEnd - mYStart);
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
        if(mCurrentY == mYEnd) {
            return 0;
        }
        if(mSampleBuffer == NULL) {
            // 2(imageX, imageY) + quota size(extra requested 1/2D samples)
            size_t bufferSize = mSamplesPerPixel * (2 + mSampleQuota.size());
            mSampleBuffer = new float[bufferSize];
        }

        float* imageBuffer = mSampleBuffer;
        float* quotaBuffer = mSampleBuffer + 2 * mSamplesPerPixel;
        // 1 pixel is 1 strata, gettting 2d stratified as 
        // mSamplesPerPixel sub strata
        stratifiedUniform2D(imageBuffer, 1);
        float* currentOffset = quotaBuffer;
        for(size_t i = 0; i < mSampleQuota.n1D.size(); ++i) {
            stratifiedUniform1D(currentOffset, mSampleQuota.n1D[i]);
            currentOffset += mSampleQuota.n1D[i] * mSamplesPerPixel;
        }
        for(size_t i = 0; i < mSampleQuota.n2D.size(); ++i) {
            stratifiedUniform2D(currentOffset, mSampleQuota.n2D[i]);
            currentOffset += 2 * mSampleQuota.n2D[i] * mSamplesPerPixel;
        }
        // shuffle the above stratified result
        float* shuffleBuffer = quotaBuffer;
        for(size_t i = 0; i < mSampleQuota.n1D.size(); ++i) {
            for(uint32_t j = 0; j < mSampleQuota.n1D[i]; ++j) {
                shuffle(shuffleBuffer, mSamplesPerPixel, 1);
                shuffleBuffer += mSamplesPerPixel;
            }
        }
        for(size_t i = 0; i < mSampleQuota.n2D.size(); ++i) {
            for(uint32_t j = 0; j < mSampleQuota.n2D[i]; ++j) {
                shuffle(shuffleBuffer, mSamplesPerPixel, 2);
                shuffleBuffer += 2 * mSamplesPerPixel;
            }
        }

        // fill in mSampleBuffer to samples
        for(int i = 0; i < mSamplesPerPixel; ++i) {
            samples[i].imageX = mCurrentX + imageBuffer[2 * i];
            samples[i].imageY = mCurrentY + imageBuffer[2 * i + 1];
        }
        float* fillinBuffer = quotaBuffer;
        for(size_t i = 0; i < mSampleQuota.n1D.size(); ++i) {
            for(uint32_t j = 0; j < mSampleQuota.n1D[i]; ++j) {
                for(int k = 0; k < mSamplesPerPixel; ++k) {
                    samples[k].u1D[i][j] = fillinBuffer[k];
                }
                fillinBuffer += mSamplesPerPixel;
            }
        }
        for(size_t i = 0; i < mSampleQuota.n2D.size(); ++i) {
            for(uint32_t j = 0; j < mSampleQuota.n2D[i]; ++j) {
                for(int k = 0; k < mSamplesPerPixel; ++k) {
                    samples[k].u2D[i][2 * j] = fillinBuffer[2 * k];
                    samples[k].u2D[i][2 * j + 1] = fillinBuffer[2 * k + 1];
                }
                fillinBuffer += 2 * mSamplesPerPixel;
            }
        }

        //debugOutput(samples);
        if(++mCurrentX == mXEnd) {
            mCurrentX= mXStart;
            mCurrentY++;
        }
        return mSamplesPerPixel;
    }

    void Sampler::debugOutput(Sample* samples) {
        static bool debugFlip = true;
        if(debugFlip) {
            debugFlip = false;
            // debug output, remove me later
            float* currentOffset = mSampleBuffer;
            std::cout << "image sample \n";
            for(int i = 0; i < mSamplesPerPixel; ++i) {
                float x = currentOffset[2 * i];
                float y = currentOffset[2 * i + 1];
                std::cout << "(" << x << ", " << y << ") ";
            }
            std::cout << "\nintegrator samples\n";
            currentOffset += 2 * mSamplesPerPixel;

            for(size_t i = 0; i < mSampleQuota.n1D.size(); ++i) {
                for(size_t j = 0; j < mSampleQuota.n1D[i]; ++j) {
                    for(int k = 0; k < mSamplesPerPixel; ++k) {
                        std::cout << currentOffset[k] << " ";
                    }
                    std::cout << std::endl;
                    currentOffset += mSamplesPerPixel;
                }
                std::cout << std::endl;
            }

            for(size_t i = 0; i < mSampleQuota.n2D.size(); ++i) {
                int root = (int)sqrtf((float)mSampleQuota.n2D[i]);
                for(size_t j = 0; j < mSampleQuota.n2D[i]; ++j) {
                    for(int k = 0; k < mSamplesPerPixel; ++k) {
                        int index = 2 * k;
                        std::cout <<"(" << currentOffset[index] << 
                            ", "<< currentOffset[index + 1] << ") ";
                    }
                    currentOffset += 2 * mSamplesPerPixel;
                    std::cout << std::endl;
                }
                std::cout << std::endl;
            }

            for(int i = 0; i < mSamplesPerPixel; ++i) {
                std::cout << "samples " << i << std::endl;
                std::cout << "n1d\n";
                for(size_t j = 0; j < samples[i].n1D.size(); ++j) {
                    for(size_t k = 0; k < samples[i].n1D[j]; ++k) {
                        std::cout << samples[i].u1D[j][k] << " ";
                    }
                    std::cout << std::endl;
                }
                std::cout << "n2d\n";
                for(size_t j = 0; j < samples[i].n2D.size(); ++j) {
                    for(size_t k = 0; k < samples[i].n2D[j]; ++k) {
                        std::cout << "(" << samples[i].u2D[j][2 * k] << 
                            ", " << samples[i].u2D[j][2 * k + 1] << ") ";
                    }
                    std::cout << std::endl;
                }
            }
        }
    }

    uint32_t Sampler::requestOneDQuota(uint32_t samplesNum) {
        mSampleQuota.n1D.push_back(roundToSquare(samplesNum));
        return mSampleQuota.n1D.size() - 1;
    }

    uint32_t Sampler::requestTwoDQuota(uint32_t samplesNum) {
        mSampleQuota.n2D.push_back(roundToSquare(samplesNum));
        return mSampleQuota.n2D.size() - 1;
    }

    Sample* Sampler::allocateSampleBuffer(size_t bufferSize) {
        Sample* samples = new Sample[bufferSize];
        for(size_t i = 0; i < bufferSize; ++i) {
            samples[i].allocateQuota(mSampleQuota);
        }
        return samples;
    }

    void Sampler::stratifiedUniform1D(float* buffer, uint32_t n1D) {
        float strataSize = 1.0f / (float)n1D;
        float subStrataSize = strataSize / mSamplesPerPixel;
        for(uint32_t i = 0; i < n1D; ++i) {
            for(int j = 0; j < mSamplesPerPixel; ++j) {
                float nOffset = mJitter ? j + randomFloat() : j + 0.5f;
                int index = i * mSamplesPerPixel + j;
                buffer[index] = i * strataSize + nOffset * subStrataSize;
            }
        }
    }

    void Sampler::stratifiedUniform2D(float* buffer, uint32_t n2D) {
        int root = (int)sqrtf((float)n2D);
        float strataSize = 1.0f / root;
        float subStrataSize = strataSize / mXPerPixel;
        for(uint32_t n = 0; n < n2D; ++n) {
            int uX = n % root;
            int uY = n / root;
            for(int p = 0; p < mSamplesPerPixel; ++p) {
                int pX = p % mXPerPixel;
                int pY = p / mXPerPixel;
                float xOffset = mJitter ? pX + randomFloat() : pX + 0.5f;
                float yOffset = mJitter ? pY + randomFloat() : pY + 0.5f;
                int index = 2 * (n * mSamplesPerPixel + pY * mXPerPixel + pX);
                buffer[index] = 
                    uX * strataSize + xOffset * subStrataSize; 
                buffer[index + 1 ] = 
                    uY * strataSize + yOffset * subStrataSize; 
            }
        }
    }

    void Sampler::shuffle(float* buffer, uint32_t num, uint32_t dim) {
        for(uint32_t n = 0; n < num; ++n) {
            size_t toShuffle = randomUInt() % num;
            for(uint32_t d = 0; d < dim; ++d) {
                swap(buffer[n * dim + d], buffer[toShuffle * dim + d]);
            }
        }
    }


    CDF1D::CDF1D(const vector<float>& f1D): mFunction(f1D) {
        size_t n = f1D.size();
        mDx = 1.0f / n;
        mCDF.resize(n + 1);
        mCDF[0] = 0.0f;
        // accumulate up the integral
        for(size_t i = 1; i < n + 1; ++i) {
            mCDF[i] = mCDF[i - 1] + (mFunction[i - 1] * mDx);
        }
        mIntegral = mCDF[n];
        // normalize CDF with integral
        for(size_t i = 1; i < n + 1; ++i) {
            mCDF[i] /= mIntegral;
        }
    }

    int CDF1D::sampleDiscrete(float u , float* pdf) {
        vector<float>::iterator lowBound;
        lowBound = std::lower_bound(mCDF.begin(), mCDF.end(), u);
        int offset = max(0, lowBound - mCDF.begin() - 1);
        if(pdf) {
            *pdf = (mFunction[offset] / mIntegral) * mDx;
        }
        return offset;
    }
}