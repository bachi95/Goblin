#include "GoblinFilm.h"
#include "GoblinFilter.h"
#include "GoblinUtils.h"
#include "GoblinSampler.h"
#include "GoblinImageIO.h"

#include <cstring>

namespace Goblin {

    ImageTile::ImageTile(int tileWidth, int rowId, int rowNum, 
        int colId, int colNum, const ImageRect& imageRect, 
        const Filter* filter, const float* filterTable,
        bool requireLightMap):
        mTileWidth(tileWidth), mRowId(rowId), mRowNum(rowNum),
        mColId(colId), mColNum(colNum), mImageRect(imageRect), mPixels(NULL),  
        mFilter(filter), mFilterTable(filterTable),
        mInvPixelArea(1.0f), mTotalSampleCount(0) {
        if (requireLightMap) {
            // when rendering involve photon splatting on film, we need to
            // allocate the full size image for each tile since photon can
            // splat on any location on the film instead of just one
            // particular square zone
            mTileRect = mImageRect;
            mPixels = new Pixel[mImageRect.xCount * mImageRect.yCount];
        } else {
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

    void ImageTile::addSample(float imageX, float imageY, const Color& L) {
        if(L.isNaN()) {
            cout << "sample ("<< imageX << " " << imageY
                << ") generate NaN point, discard this sample" << endl;
            return;
        }
        // transform continuous space sample to discrete space
        float dImageX = imageX - 0.5f;
        float dImageY = imageY - 0.5f;
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
                float weight = mFilterTable[iy * FILTER_TABLE_WIDTH + ix] *
                    mInvPixelArea;
                int index = (y - mTileRect.yStart) * mTileRect.xCount + 
                    (x - mTileRect.xStart);
                mPixels[index].color += weight * L;
                mPixels[index].weight += weight; 
            }
        }
    }

    void ImageTile::setTotalSampleCount(uint64_t totalSampleCount) {
        mTotalSampleCount = totalSampleCount;
    }

    uint64_t ImageTile::getTotalSampleCount() const {
        return mTotalSampleCount;
    }

    Film::Film(int xRes, int yRes, const float crop[4],
        Filter* filter, const std::string& filename,
        bool requireLightMap, bool toneMapping,
        float bloomRadius, float bloomWeight):
        mXRes(xRes), mYRes(yRes), mFilter(filter), mFilename(filename),
        mTileWidth(64), mToneMapping(toneMapping),
        mBloomRadius(bloomRadius), mBloomWeight(bloomWeight),
        mInvPixelArea(1.0f) {

        memcpy(mCrop, crop, 4 * sizeof(float));

        mXStart = ceilInt(mXRes * mCrop[0]);
        mXCount = max(1, ceilInt(mXRes * mCrop[1]) - mXStart);
        mYStart = ceilInt(mYRes * mCrop[2]); 
        mYCount = max(1, ceilInt(mYRes * mCrop[3]) - mYStart);

        mPixels = new Pixel[mXRes * mYRes];

        mInvXRes = 1.0f / (float)mXRes;
        mInvYRes = 1.0f / (float)mYRes;
        // precompute filter equation as a lookup table
        // we only computer the right up corner since for all the supported
        // filters: f(x, y) = f(|x|, |y|)
        float deltaX = mFilter->getXWidth() / FILTER_TABLE_WIDTH;
        float deltaY = mFilter->getYWidth() / FILTER_TABLE_WIDTH;
        float normalizeTerm = mFilter->getNormalizeTerm();
        size_t index = 0;
        for(int y = 0; y < FILTER_TABLE_WIDTH; ++y) {
            float fy = y * deltaY;
            for(int x = 0; x < FILTER_TABLE_WIDTH; ++x) {
                float fx = x * deltaX;
                mFilterTable[index++] = mFilter->evaluate(fx, fy) /
                    normalizeTerm;
            }
        }

        // construc image tiles
        ImageRect r(mXStart, mYStart, mXCount, mYCount);
        int colNum = ceilInt((float)mXCount / (float)mTileWidth);
        int rowNum = ceilInt((float)mYCount / (float)mTileWidth);
        for(int i = 0; i < rowNum; ++i) {
            for(int j = 0; j < colNum; ++j) {
                ImageTile* tile = new ImageTile(mTileWidth, 
                    i, rowNum, j, colNum, r, mFilter, mFilterTable,
                    requireLightMap);
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

    void Film::addSample(float imageX, float imageY, const Color& L) {
        if(L.isNaN()) {
            cout << "sample ("<< imageX << " " << imageY
                << ") generate NaN point, discard this sample" << endl;
            return;
        }
        // transform continuous space sample to discrete space
        float dImageX = imageX - 0.5f;
        float dImageY = imageY - 0.5f;
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
                float weight = mFilterTable[iy * FILTER_TABLE_WIDTH + ix] *
                    mInvPixelArea;
                int index = y * mXRes + x;
                mPixels[index].color += weight * L;
                mPixels[index].weight += weight; 
            }
        }
    } 

    void Film::setFilmArea(float filmArea) {
        mFilmArea = filmArea;
        mInvPixelArea = mXRes * mYRes / filmArea;
        for (size_t i = 0; i < mTiles.size(); ++i) {
            mTiles[i]->setInvPixelArea(mInvPixelArea);
        }
    }

    void Film::mergeTiles() {
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
    }

    void Film::scaleImage(float scale) {
        for (int y = 0; y < mYRes; ++y) {
            for (int x = 0; x < mXRes; ++x) {
                int index = mXRes * y + x;
                mPixels[index].color *= scale;
            }
        }
    }

    void Film::writeImage(bool normalize) {
        Color* colors = new Color[mXRes * mYRes];
        for(int y = 0; y < mYRes; ++y) {
            for(int x = 0; x < mXRes; ++x) {
                int index = mXRes * y + x;
                colors[index] = normalize?
                    mPixels[index].color / mPixels[index].weight:
                    mPixels[index].color;
            }
        }

        // draw debug info
        cout << "drawing debug info: " << endl;
        for(size_t i = 0; i < mTiles.size(); ++i) {
            const DebugInfo& debugInfo = mTiles[i]->getDebugInfo();
            const vector<pair<DebugLine, Color> >& lines = 
                debugInfo.getLines();
            for(size_t i = 0; i < lines.size(); ++i) {
                const DebugLine& line = lines[i].first;
                drawLine(line.first, line.second, colors, 
                    mXRes, mYRes, lines[i].second);
            }
            const vector<pair<Vector2, Color> >& points = 
                debugInfo.getPoints();
            for(size_t i = 0; i < points.size(); ++i) {
                drawPoint(points[i].first, colors, mXRes, mYRes, 
                    points[i].second, 1);
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
        bool requireLightMap = params.getBool("require_light_map");
        return new Film(xRes, yRes, crop, filter, filePath, 
            requireLightMap, toneMapping, bloomRadius, bloomWeight);
    }

}

