#include "GoblinFilm.h"
#include "GoblinFilter.h"
#include "GoblinUtils.h"
#include "GoblinSampler.h"
#include "GoblinImageIO.h"

#include <cstring>

namespace Goblin {

    ImageTile::ImageTile(int tileWidth, int rowId, int rowNum, 
        int colId, int colNum, const ImageRect& imageRect, 
        const Filter* filter, const float* filterTable):
        mTileWidth(tileWidth), mRowId(rowId), mRowNum(rowNum),
        mColId(colId), mColNum(colNum), mImageRect(imageRect), mPixels(NULL),  
        mFilter(filter), mFilterTable(filterTable) {
        int xStart = mImageRect.xStart + mTileWidth * mColId;
        int yStart = mImageRect.yStart + mTileWidth * mRowId;
        int xEnd = xStart + mTileWidth;
        int yEnd = yStart + mTileWidth;
        if(mColId == mColNum - 1) {
            xEnd = min(mImageRect.xStart + mImageRect.xCount, xEnd);
        }
        if(mRowId == mRowNum - 1) {
            yEnd = min(mImageRect.yStart + mImageRect.yCount, yEnd);
        }
        // neighbor tiles will share an area that both get updated
        // because of the image filter, make each one has its own
        // local area of that share part so we don't need to worry
        // about concurrency issue
        float filterXWidth = mFilter->getXWidth();
        float filterYWidth = mFilter->getYWidth();
        if(mColId != 0) {
            xStart = floorInt(xStart + 0.5f - filterXWidth);
        }
        if(mRowId != 0) {
            yStart = floorInt(yStart + 0.5f - filterYWidth);
        }
        if(mColId != colNum - 1) {
            xEnd = floorInt(xEnd + 0.5f + filterXWidth);
        }
        if(mRowId != rowNum - 1) {
            yEnd = floorInt(yEnd + 0.5f + filterYWidth);
        }

        int xCount = xEnd - xStart;
        int yCount = yEnd - yStart;
        mTileRect.xStart = xStart;
        mTileRect.yStart = yStart;
        mTileRect.xCount = xCount;
        mTileRect.yCount = yCount;
        mPixels = new Pixel[xCount * yCount];
    }

    ImageTile::~ImageTile() {
        if(mPixels) {
            delete [] mPixels;
            mPixels = NULL;
        }
    }

    void ImageTile::getImageRange(int* xStart, int *xEnd,
        int* yStart, int *yEnd) const {
        *xStart = mTileRect.xStart;
        *xEnd = mTileRect.xStart + mTileRect.xCount;
        *yStart = mTileRect.yStart;
        *yEnd = mTileRect.yStart + mTileRect.yCount;
    }

    void ImageTile::getSampleRange(int* xStart, int* xEnd, 
        int* yStart, int* yEnd) const {
        float xWidth = mFilter->getXWidth();
        float yWidth = mFilter->getYWidth();

        *xStart = mImageRect.xStart + mTileWidth * mColId;
        *yStart = mImageRect.yStart + mTileWidth * mRowId;
        *xEnd = *xStart + mTileWidth;
        *yEnd = *yStart + mTileWidth;
        if(mColId == mColNum - 1) {
            *xEnd = min(mImageRect.xStart + mImageRect.xCount, *xEnd);
        }
        if(mRowId == mRowNum - 1) {
            *yEnd = min(mImageRect.yStart + mImageRect.yCount, *yEnd);
        }
 
        // border case need to take some sample outside of the image
        if(mColId == 0) {
            *xStart = floorInt(*xStart + 0.5f - xWidth);
        }
        if(mColId == mColNum - 1) {
            *xEnd = floorInt(*xEnd + 0.5f + xWidth);
        }
        if(mRowId == 0) {
            *yStart = floorInt(*yStart + 0.5f - yWidth);
        } 
        if(mRowId == mRowNum - 1) {
            *yEnd = floorInt(*yEnd + 0.5f + yWidth);
        }
    }

    const Pixel* ImageTile::getTileBuffer() const {
        return mPixels;
    }

    void ImageTile::addSample(const Sample& sample, const Color& L) {
        if(L.isNaN()) {
            cout << "sample ("<< sample.imageX << " " << sample.imageY
                << ") generate NaN point, discard this sample" << endl;
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
        x0 = max(x0, mTileRect.xStart);
        x1 = min(x1, mTileRect.xStart + mTileRect.xCount - 1);
        y0 = max(y0, mTileRect.yStart);
        y1 = min(y1, mTileRect.yStart + mTileRect.yCount - 1);
        
        for(int y = y0; y <= y1; ++y) {
            float fy = fabs(FILTER_TABLE_WIDTH * (y - dImageY) / yWidth);
            int iy = min(floorInt(fy), FILTER_TABLE_WIDTH - 1);
            for(int x = x0; x <= x1; ++x) {
                float fx = fabs(FILTER_TABLE_WIDTH * (x - dImageX) / xWidth);
                int ix = min(floorInt(fx), FILTER_TABLE_WIDTH - 1);
                float weight = mFilterTable[iy * FILTER_TABLE_WIDTH + ix];
                int index = (y - mTileRect.yStart) * mTileRect.xCount + 
                    (x - mTileRect.xStart);
                mPixels[index].color += weight * L;
                mPixels[index].weight += weight; 
            }
        }
    } 


    Film::Film(int xRes, int yRes, const float crop[4],
        Filter* filter, const std::string& filename,
        int tileWidth, bool toneMapping, 
        float bloomRadius, float bloomWeight): 
        mXRes(xRes), mYRes(yRes), mFilter(filter), mFilename(filename),
        mTileWidth(tileWidth), mToneMapping(toneMapping), 
        mBloomRadius(bloomRadius), mBloomWeight(bloomWeight) {

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

        // construc image tiles
        ImageRect r(mXStart, mYStart, mXCount, mYCount);
        int colNum = ceilInt((float)mXCount / (float)mTileWidth);
        int rowNum = ceilInt((float)mYCount / (float)mTileWidth);
        for(int i = 0; i < rowNum; ++i) {
            for(int j = 0; j < colNum; ++j) {
                ImageTile* tile = new ImageTile(mTileWidth, 
                    i, rowNum, j, colNum, r, mFilter, mFilterTable);
                mTiles.push_back(tile);
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

        for(size_t i = 0; i < mTiles.size(); ++i) {
            delete mTiles[i];
        } 
        mTiles.clear();
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
            cout << "sample ("<< sample.imageX << " " << sample.imageY
                << ") generate NaN point, discard this sample" << endl;
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
        // merget tiles
        for(size_t i = 0; i < mTiles.size(); ++i) {
            int xStart, xEnd, yStart, yEnd;
            mTiles[i]->getImageRange(&xStart, &xEnd, &yStart, &yEnd);
            const Pixel* tileBuffer = mTiles[i]->getTileBuffer();
            int tileWidth = xEnd - xStart;
            for(int y = yStart; y < yEnd; ++y) {
                for(int x = xStart; x < xEnd; ++x) {
                    int tileIndex = (y - yStart) * tileWidth + (x - xStart);
                    int filmIndex = y * mXRes + x;
                    mPixels[filmIndex].color += tileBuffer[tileIndex].color;
                    mPixels[filmIndex].weight += tileBuffer[tileIndex].weight;
                }
            }
        }

        Color* colors = new Color[mXRes * mYRes];
        for(int y = 0; y < mYRes; ++y) {
            for(int x = 0; x < mXRes; ++x) {
                int index = mXRes * y + x;
                colors[index] = mPixels[index].color / mPixels[index].weight;
            }
        }
        cout << "write image to : " << mFilename << endl;
        if(mBloomRadius > 0.0f && mBloomWeight > 0.0f) {
            Goblin::bloom(colors, mXRes, mYRes, mBloomRadius, mBloomWeight);
        }
        Goblin::writeImage(mFilename, colors, mXRes, mYRes, mToneMapping);
        delete [] colors;
    }

    Film* ImageFilmCreator::create(const ParamSet& params, 
        Filter* filter) const {
        Vector2 res = params.getVector2("resolution", Vector2(640, 480));
        int xRes = static_cast<int>(res.x);
        int yRes = static_cast<int>(res.y);

        Vector4 windowCrop = params.getVector4("crop", Vector4(0, 1, 0, 1));
        float crop[4];
        for(int i = 0; i < 4; ++i) {
            crop[i] = windowCrop[i];
        }
        string filePath = params.getString("file", "goblin.png");
        bool toneMapping = params.getBool("tone_mapping");
        float bloomRadius = params.getFloat("bloom_radius");
        float bloomWeight = params.getFloat("bloom_weight");
        int tileWidth = params.getInt("tile_width", 16);
        return new Film(xRes, yRes, crop, filter, filePath, 
            tileWidth, toneMapping, bloomRadius, bloomWeight);
    }

}

