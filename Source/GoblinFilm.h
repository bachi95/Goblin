#ifndef GOBLIN_FILM_H
#define GOBLIN_FILM_H

#include "GoblinColor.h"
#include "GoblinFactory.h"
#include "GoblinUtils.h"
#include "GoblinVector.h"

namespace Goblin {

    const int FILTER_TABLE_WIDTH = 16;
    class Sample;
    class Filter;

    class Pixel {
    public:
        Pixel(): color(Color::Black), weight(0.0f) {}
        Color color;
        float weight;
        float pad[3];
    };

    typedef pair<Vector2, Vector2> DebugLine;

    class DebugInfo {
    public:
        DebugInfo() {}
        void addLine(const DebugLine& line);
        void addPoint(const Vector2& point);
        const vector<DebugLine>& getLines() const;
        const vector<Vector2>& getPoints() const;
    private:
        vector<DebugLine> mLines;
        vector<Vector2> mPoints;
    };

    inline void DebugInfo::addLine(const DebugLine& line) { 
        mLines.push_back(line); 
    }

    inline void DebugInfo::addPoint(const Vector2& point) { 
        mPoints.push_back(point); 
    }

    inline const vector<DebugLine>& DebugInfo::getLines() const {
        return mLines;
    }

    inline const vector<Vector2>& DebugInfo::getPoints() const {
        return mPoints;
    }

    struct ImageRect {
        ImageRect() {}
        ImageRect(int x, int y, int w, int h): 
            xStart(x), yStart(y), xCount(w), yCount(h) {}
        int xStart, yStart, xCount, yCount;
    };

    class ImageTile {
    public:
        ImageTile(int tileWidth, int rowId, int rowNum, int colId, int colNum, 
            const ImageRect& imageRect, const Filter* filter, 
            const float* filterTable);
        ~ImageTile();
        void getImageRange(int* xStart, int *xEnd,
            int* yStart, int* yEnd) const;
        void getSampleRange(int* xStart, int* xEnd,
            int* yStart, int* yEnd) const;
        const Pixel* getTileBuffer() const;
        const DebugInfo& getDebugInfo() const;
        void addSample(const Sample& sample, const Color& L);
        void addDebugLine(const DebugLine& line);
        void addDebugPoint(const Vector2& point);
    private:
        int mTileWidth;
        int mRowId, mRowNum;
        int mColId, mColNum;
        ImageRect mTileRect;
        ImageRect mImageRect;
        Pixel* mPixels;
        const Filter* mFilter;
        const float* mFilterTable;
        DebugInfo mDebugInfo;
    };

    inline const Pixel* ImageTile::getTileBuffer() const {
        return mPixels;
    }

    inline const DebugInfo& ImageTile::getDebugInfo() const {
        return mDebugInfo;
    }

    inline void ImageTile::addDebugLine(const DebugLine& line) {
        mDebugInfo.addLine(line);
    }

    inline void ImageTile::addDebugPoint(const Vector2& point) {
        mDebugInfo.addPoint(point);
    }


    class Film {
    public:
        Film(int xRes, int yRes, const float crop[4], 
            Filter* filter, const std::string& filename, 
            int tileWidth, bool toneMapping = false, 
            float bloomRadius = 0.0f, float bloomWeight = 0.0f);
        ~Film();

        int getXResolution() const;
        int getYResolution() const;
        void getSampleRange(int* xStart, int* xEnd,
            int* yStart, int* yEnd) const;
        vector<ImageTile*>& getTiles();
        void addSample(const Sample& sample, const Color& L);
        void writeImage();

    private:
        int mXRes, mYRes;
        int mXStart, mYStart, mXCount, mYCount;
        float mFilterTable[FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH];
        float mCrop[4];
        Filter* mFilter;
        Pixel* mPixels;
        vector<ImageTile*> mTiles;
        std::string mFilename;
        int mTileWidth;
        bool mToneMapping;
        float mBloomRadius;
        float mBloomWeight;
    };

    inline int Film::getXResolution() const { return mXRes; }
    inline int Film::getYResolution() const { return mYRes; }
    inline vector<ImageTile*>& Film::getTiles() { return mTiles; }

    class ImageFilmCreator : public Creator<Film, const ParamSet&, Filter*> {
    public:
        Film* create(const ParamSet& params, Filter* filter) const;
    };
}

#endif //GOBLIN_FILM_H
