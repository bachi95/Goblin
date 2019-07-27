#ifndef GOBLIN_FILM_H
#define GOBLIN_FILM_H

#include "GoblinColor.h"
#include "GoblinDebugData.h"
#include "GoblinUtils.h"
#include "GoblinVector.h"

namespace Goblin {

const int FILTER_TABLE_WIDTH = 16;
class Sample;
struct SampleRange;
class Filter;

class Pixel {
public:
    Pixel(): color(Color::Black), weight(0.0f) {}
    Color color;
    float weight;
    float pad[3];
};

struct ImageRect {
    ImageRect() {}
    ImageRect(int x, int y, int w, int h):
        xStart(x), yStart(y), xCount(w), yCount(h) {}

    int pixelNum() const { return xCount * yCount; }

    int pixelToOffset(int x, int y) const {
        return  (y - yStart) * xCount + (x - xStart);
    }

    void offsetToPixel(int offset, int* x, int* y) const {
        *x = offset % xCount + xStart;
        *y = offset / xCount + yStart;
    }

    int xStart, yStart, xCount, yCount;
};

class FilterTable {
public:
    FilterTable(const Filter* filter);

    float evaluate(float x, float y) const;

    const Vector2& getFilterWidth() const {
        return mFilterWidth;
    }

private:
    float mTable[FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH];
    Vector2 mFilterWidth;
};

class ImageTile {
public:
    ImageTile(const ImageRect& tileRect, const FilterTable& cachedFilter);

    ~ImageTile();

    void getTileRange(int* xStart, int *xEnd,
        int* yStart, int* yEnd) const;

    const Pixel* getTileBuffer() const;

    void addSample(float imageX, float imageY, const Color& L);

private:
    ImageRect mTileRect;
    Pixel* mPixels;
    const FilterTable& mCachedFilter;
};

inline const Pixel* ImageTile::getTileBuffer() const {
    return mPixels;
}

class Film {
public:
    Film(int xRes, int yRes, const float crop[4],
        Filter* filter, const std::string& filename,
        bool toneMapping = false,
        float bloomRadius = 0.0f, float bloomWeight = 0.0f);

    ~Film();

    int getXResolution() const;

    int getYResolution() const;

    float getInvXResolution() const;

    float getInvYResolution() const;

    void getImageRect(ImageRect& imageRect) const;

    void getSampleRange(SampleRange& sampleRange) const;

    const FilterTable& getFilterTable() const {
        return mCachedFilter;
    }

    void scaleImage(float s);

    void writeImage(bool normalize = true);

    void mergeTile(const ImageTile& tile);

    void addDebugLine(const DebugLine& l, const Color& c);

    void addDebugPoint(const Vector2& p, const Color& c);

private:
    int mXRes, mYRes;
    int mXStart, mYStart, mXCount, mYCount;
    float mInvXRes, mInvYRes;
    float mCrop[4];
    Filter* mFilter;
    FilterTable mCachedFilter;
    Pixel* mPixels;
    std::string mFilename;
    bool mToneMapping;
    float mBloomRadius;
    float mBloomWeight;
    std::vector<std::pair<DebugLine, Color> > mDebugLines;
    std::vector<std::pair<Vector2, Color> > mDebugPoints;
};

inline int Film::getXResolution() const { return mXRes; }

inline int Film::getYResolution() const { return mYRes; }

inline float Film::getInvXResolution() const { return mInvXRes; }
    
inline float Film::getInvYResolution() const { return mInvYRes; }

Film* createImageFilm(const ParamSet& params, Filter* filter);

}

#endif //GOBLIN_FILM_H
