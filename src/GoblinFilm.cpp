#include "GoblinFilm.h"
#include "GoblinFilter.h"
#include "GoblinUtils.h"
#include "GoblinSampler.h"
#include "GoblinImageIO.h"

#include <cstring>

namespace Goblin {
FilterTable::FilterTable(const Filter* filter):
    mFilterWidth(filter->getXWidth(), filter->getYWidth()) {
    // precompute filter equation as a lookup table
    // we only computer the right up corner since for all the supported
    // filters: f(x, y) = f(|x|, |y|)
    float deltaX = mFilterWidth.x / FILTER_TABLE_WIDTH;
    float deltaY = mFilterWidth.y / FILTER_TABLE_WIDTH;
    float normalizeTerm = filter->getNormalizeTerm();
    size_t index = 0;
    for (int y = 0; y < FILTER_TABLE_WIDTH; ++y) {
        float fy = y * deltaY;
        for (int x = 0; x < FILTER_TABLE_WIDTH; ++x) {
            float fx = x * deltaX;
            mTable[index++] = filter->evaluate(fx, fy) /
                normalizeTerm;
        }
    }
}
     
float FilterTable::evaluate(float x, float y) const {
    int iy = std::min(
        floorInt(fabs(FILTER_TABLE_WIDTH * y / mFilterWidth.y)),
        FILTER_TABLE_WIDTH - 1);
    int ix = std::min(
        floorInt(fabs(FILTER_TABLE_WIDTH * x / mFilterWidth.x)),
        FILTER_TABLE_WIDTH - 1);
    return mTable[iy * FILTER_TABLE_WIDTH + ix];
}

ImageTile::ImageTile(const ImageRect& tileRect,
    const FilterTable& cachedFilter):
    mTileRect(tileRect), mPixels(nullptr),
    mCachedFilter(cachedFilter) {
    mPixels = new Pixel[mTileRect.xCount * mTileRect.yCount];
}

ImageTile::~ImageTile() {
    if (mPixels) {
        delete [] mPixels;
        mPixels = nullptr;
    }
}

void ImageTile::getTileRange(int* xStart, int *xEnd,
    int* yStart, int *yEnd) const {
    *xStart = mTileRect.xStart;
    *xEnd = mTileRect.xStart + mTileRect.xCount;
    *yStart = mTileRect.yStart;
    *yEnd = mTileRect.yStart + mTileRect.yCount;
}

void ImageTile::addSample(float imageX, float imageY, const Color& L) {
    if (L.isNaN()) {
		std::cout << "sample ("<< imageX << " " << imageY
            << ") generate NaN point, discard this sample" << std::endl;
        return;
    }
    // transform continuous space sample to discrete space
    float dImageX = imageX - 0.5f;
    float dImageY = imageY - 0.5f;
    // calculate the pixel range covered by filter center at sample
    const Vector2 filterWidth = mCachedFilter.getFilterWidth();
    int x0 = ceilInt(dImageX - filterWidth.x);
    int x1 = floorInt(dImageX + filterWidth.x);
    int y0 = ceilInt(dImageY - filterWidth.y);
    int y1 = floorInt(dImageY + filterWidth.y);
    x0 = std::max(x0, mTileRect.xStart);
    x1 = std::min(x1, mTileRect.xStart + mTileRect.xCount - 1);
    y0 = std::max(y0, mTileRect.yStart);
    y1 = std::min(y1, mTileRect.yStart + mTileRect.yCount - 1);

    for (int y = y0; y <= y1; ++y) {
        for (int x = x0; x <= x1; ++x) {
            float w = mCachedFilter.evaluate(x - dImageX, y - dImageY);
            int index = (y - mTileRect.yStart) * mTileRect.xCount +
                (x - mTileRect.xStart);
            mPixels[index].color += w * L;
            mPixels[index].weight += w;
        }
    }
}

Film::Film(int xRes, int yRes, const float crop[4],
    Filter* filter, const std::string& filename,
    bool toneMapping,
    float bloomRadius, float bloomWeight):
    mXRes(xRes), mYRes(yRes), mFilter(filter), mCachedFilter(filter),
    mFilename(filename), mToneMapping(toneMapping),
    mBloomRadius(bloomRadius), mBloomWeight(bloomWeight) {

    memcpy(mCrop, crop, 4 * sizeof(float));

    mXStart = ceilInt(mXRes * mCrop[0]);
    mXCount = std::max(1, ceilInt(mXRes * mCrop[1]) - mXStart);
    mYStart = ceilInt(mYRes * mCrop[2]);
    mYCount = std::max(1, ceilInt(mYRes * mCrop[3]) - mYStart);

    mPixels = new Pixel[mXRes * mYRes];

    mInvXRes = 1.0f / (float)mXRes;
    mInvYRes = 1.0f / (float)mYRes;
}
    
Film::~Film() {
    if (mPixels != nullptr) {
        delete[] mPixels;
        mPixels = nullptr;
    }
    if (mFilter != nullptr) {
        delete mFilter;
        mFilter = nullptr;
    }
}

void Film::getImageRect(ImageRect& imageRect) const {
    imageRect.xStart = mXStart;
    imageRect.yStart = mYStart;
    imageRect.xCount = mXCount;
    imageRect.yCount = mYCount;
}

void Film::getSampleRange(SampleRange& sampleRange) const {
    float xWidth = mFilter->getXWidth();
    float yWidth = mFilter->getYWidth();
    sampleRange.xStart = floorInt(mXStart + 0.5f - xWidth);
    sampleRange.xEnd = floorInt(mXStart + 0.5f + mXCount + xWidth);
    sampleRange.yStart = floorInt(mYStart + 0.5f - yWidth);
    sampleRange.yEnd = floorInt(mYStart + 0.5f + mYCount + yWidth);
}

void Film::mergeTile(const ImageTile& tile) {
    int xStart, xEnd, yStart, yEnd;
    tile.getTileRange(&xStart, &xEnd, &yStart, &yEnd);
    const Pixel* tileBuffer = tile.getTileBuffer();
    int tileWidth = xEnd - xStart;
    for (int y = yStart; y < yEnd; ++y) {
        for (int x = xStart; x < xEnd; ++x) {
            int tileIndex = (y - yStart) * tileWidth + (x - xStart);
            int filmIndex = y * mXRes + x;
            mPixels[filmIndex].color += tileBuffer[tileIndex].color;
            mPixels[filmIndex].weight += tileBuffer[tileIndex].weight;
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
    std::vector<Color> colors(mXRes * mYRes);
    for (int y = 0; y < mYRes; ++y) {
        for (int x = 0; x < mXRes; ++x) {
            int index = mXRes * y + x;
            colors[index] = normalize?
                mPixels[index].color / mPixels[index].weight:
                mPixels[index].color;
        }
    }

    // draw debug info
	std::cout << "drawing debug info: " << std::endl;
    for (size_t i = 0; i < mDebugLines.size(); ++i) {
        const DebugLine& line = mDebugLines[i].first;
        drawLine(line.first, line.second, colors.data(),
            mXRes, mYRes, mDebugLines[i].second);
    }
    for (size_t i = 0; i < mDebugPoints.size(); ++i) {
        drawPoint(mDebugPoints[i].first, colors.data(), mXRes, mYRes,
            mDebugPoints[i].second, 1);
    }

	std::cout << "write image to : " << mFilename << std::endl;
    if (mBloomRadius > 0.0f && mBloomWeight > 0.0f) {
        Goblin::bloom(colors.data(), mXRes, mYRes, mBloomRadius, mBloomWeight);
    }
    Goblin::writeImage(mFilename, colors.data(), mXRes, mYRes, mToneMapping);
}

void Film::addDebugLine(const DebugLine& l, const Color& c) {
    mDebugLines.push_back(std::pair<DebugLine, Color>(l, c));
}

void Film::addDebugPoint(const Vector2& p, const Color& c) {
    mDebugPoints.push_back(std::pair<Vector2, Color>(p, c));
}

Film* createImageFilm(const ParamSet& params, Filter* filter) {
	Vector2 res = params.getVector2("resolution", Vector2(512, 512));
	int xRes = static_cast<int>(res.x);
	int yRes = static_cast<int>(res.y);

	Vector4 windowCrop = params.getVector4("crop", Vector4(0, 1, 0, 1));
	float crop[4];
	for (int i = 0; i < 4; ++i) {
		crop[i] = windowCrop[i];
	}
	std::string filePath = params.getString("file", "goblin.exr");
	bool toneMapping = params.getBool("tone_mapping");
	float bloomRadius = params.getFloat("bloom_radius");
	float bloomWeight = params.getFloat("bloom_weight");
	return new Film(xRes, yRes, crop, filter, filePath,
		toneMapping, bloomRadius, bloomWeight);
}

} // namespace Goblin