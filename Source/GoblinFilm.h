#ifndef GOBLIN_FILM_H
#define GOBLIN_FILM_H

#include "GoblinColor.h"
#include "GoblinFactory.h"
#include <string>
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
        void addSample(const Sample& sample, const Color& L);
    private:
        int mTileWidth;
        int mRowId, mRowNum;
        int mColId, mColNum;
        ImageRect mTileRect;
        ImageRect mImageRect;
        Pixel* mPixels;
        const Filter* mFilter;
        const float* mFilterTable;
    };

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
