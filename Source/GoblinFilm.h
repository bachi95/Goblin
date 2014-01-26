#ifndef GOBLIN_FILM_H
#define GOBLIN_FILM_H

#include "GoblinColor.h"
#include <string>
namespace Goblin {
    const int FILTER_TABLE_WIDTH = 16;
    class Sample;
    class Filter;
    class Film {
    public:
        Film(int xRes, int yRes, const float crop[4], 
            Filter* filter, const std::string& filename);
        ~Film(); 

        int getXResolution() const;
        int getYResolution() const;
        void getSampleRange(int* xStart, int* xEnd,
            int* yStart, int* yEnd) const;

        void addSample(const Sample& sample, const Color& L);
        void writeImage();

    private:
        int mXRes, mYRes;
        int mXStart, mYStart, mXCount, mYCount;
        float mFilterTable[FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH];
        float mCrop[4];
        std::string mFilename;
        class Pixel {
        public:
            Pixel(): color(Color::Black), weight(0.0f) {}
            Color color;
            float weight;
            float pad[3];
        };
        Filter* mFilter;
        Pixel* mPixels;
    };

    inline int Film::getXResolution() const { return mXRes; }
    inline int Film::getYResolution() const { return mYRes; }
}

#endif //GOBLIN_FILM_H
