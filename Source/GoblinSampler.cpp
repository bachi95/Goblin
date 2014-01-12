#include "GoblinSampler.h"

namespace Goblin {
    Sampler::Sampler(int xStart, int xEnd, int yStart, int yEnd,
        int xPerPixel, int yPerPixel):
        mXStart(xStart), mXEnd(xEnd), 
        mYStart(yStart), mYEnd(yEnd),
        mCurrentX(xStart), mCurrentY(yStart),
        mXPerPixel(xPerPixel), mYPerPixel(yPerPixel),
        mSamplesPerPixel(xPerPixel * yPerPixel) {}

    int Sampler::maxSamplesPerRequest() const {
        return mSamplesPerPixel;
    }

    int Sampler::maxTotalSamples() const {
        return mSamplesPerPixel * (mXEnd - mXStart) * (mYEnd - mYStart);
    }

    int Sampler::requestSamples(CameraSample* samples) {
        if(mCurrentY == mYEnd) {
            return 0;
        }
        stratifiedSample2D(samples);
        for(int i = 0; i < mSamplesPerPixel; ++i) {
            samples[i].imageX += mCurrentX;
            samples[i].imageY += mCurrentY;
        }
        if(++mCurrentX == mXEnd) {
            mCurrentX= mXStart;
            mCurrentY++;
        }
        return mSamplesPerPixel;
    }

    void Sampler::stratifiedSample2D(CameraSample* samples) {
        float dx = 1.0f / static_cast<float>(mXPerPixel);
        float dy = 1.0f / static_cast<float>(mYPerPixel);
        for(int y = 0; y < mYPerPixel; ++y) {
            for(int x = 0; x < mXPerPixel; ++x) {
                int index = y * mXPerPixel + x;
                samples[index].imageX = (static_cast<float>(x) + 0.5f) * dx;
                samples[index].imageY = (static_cast<float>(y) + 0.5f) * dy;
            }
        }
    }

}