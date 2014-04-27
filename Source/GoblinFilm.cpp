#include "GoblinFilm.h"
#include "GoblinFilter.h"
#include "GoblinUtils.h"
#include "GoblinSampler.h"
#include "GoblinImageIO.h"

#include <cstring>

namespace Goblin {
    Film::Film(int xRes, int yRes, const float crop[4],
        Filter* filter, const std::string& filename): 
        mXRes(xRes), mYRes(yRes), mFilter(filter), mFilename(filename) {

        memcpy(mCrop, crop, 4 * sizeof(float));

        mXStart = ceilInt(mXRes * mCrop[0]);
        mXCount = max(1, ceilInt(mXRes * mCrop[1]) - mXStart);
        mYStart = ceilInt(mYRes * mCrop[2]); 
        mYCount = max(1, ceilInt(mYRes * mCrop[3]) - mYStart);

        mPixels = new Pixel[mXRes * mYRes];

        // precompute filter equation as a lookup table
        // we only computer the right up corner since for all the supported
        // filters: f(x, y) = f(|x|, |y|)
        float deltaX = mFilter->getXWidth() / FILTER_TABLE_WIDTH;
        float deltaY = mFilter->getYWidth() / FILTER_TABLE_WIDTH;
        size_t index = 0;
        for(int y = 0; y < FILTER_TABLE_WIDTH; ++y) {
            float fy = y * deltaY;
            for(int x = 0; x < FILTER_TABLE_WIDTH; ++x) {
                float fx = x * deltaX;
                mFilterTable[index++] = mFilter->evaluate(fx, fy);
            }
        }
    } 
    
    Film::~Film() {
        if(mPixels != NULL) {
            delete[] mPixels;
            mPixels = NULL;
        }
        if(mFilter != NULL) {
            delete mFilter;
            mFilter = NULL;
        }
    }

    void Film::getSampleRange(int* xStart, int* xEnd,
        int* yStart, int* yEnd) const {
        float xWidth = mFilter->getXWidth();
        float yWidth = mFilter->getYWidth();
        *xStart = floorInt(mXStart + 0.5f - xWidth);
        *xEnd = floorInt(mXStart + 0.5f + mXCount + xWidth);
        *yStart = floorInt(mYStart + 0.5f - yWidth);
        *yEnd = floorInt(mYStart + 0.5f + mYCount + yWidth);
    }

    void Film::addSample(const Sample& sample, const Color& L) {
        if(L.isNaN()) {
            std::cout << "sample ("<< sample.imageX << " " << sample.imageY
                << ") generate NaN point, discard this sample" << std::endl;
            return;
        }
        // transform continuous space sample to discrete space
        float dImageX = sample.imageX - 0.5f;
        float dImageY = sample.imageY - 0.5f;
        // calculate the pixel range covered by filter center at sample
        float xWidth = mFilter->getXWidth();
        float yWidth = mFilter->getYWidth();
        int x0 = ceilInt(dImageX - xWidth);
        int x1 = floorInt(dImageX + xWidth);
        int y0 = ceilInt(dImageY - yWidth);
        int y1 = floorInt(dImageY + yWidth);
        x0 = max(x0, mXStart);
        x1 = min(x1, mXStart + mXCount - 1);
        y0 = max(y0, mYStart);
        y1 = min(y1, mYStart + mYCount - 1);
        
        for(int y = y0; y <= y1; ++y) {
            float fy = fabs(FILTER_TABLE_WIDTH * (y - dImageY) / yWidth);
            int iy = min(floorInt(fy), FILTER_TABLE_WIDTH - 1);
            for(int x = x0; x <= x1; ++x) {
                float fx = fabs(FILTER_TABLE_WIDTH * (x - dImageX) / xWidth);
                int ix = min(floorInt(fx), FILTER_TABLE_WIDTH - 1);
                float weight = mFilterTable[iy * FILTER_TABLE_WIDTH + ix];
                int index = y * mXRes + x;
                // TODO atomic add when this come to multi thread
                mPixels[index].color += weight * L;
                mPixels[index].weight += weight; 
            }
        }
    } 

    void Film::writeImage() {
        Color* colors = new Color[mXRes * mYRes];
        for(int y = 0; y < mYRes; ++y) {
            for(int x = 0; x < mXRes; ++x) {
                int index = mXRes * y + x;
                colors[index] = mPixels[index].color / mPixels[index].weight;
            }
        }
        std::cout << "write image to : " << mFilename << std::endl;
        Goblin::writeImage(mFilename, colors, mXRes, mYRes);
        delete [] colors;
    }

}

