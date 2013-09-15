#ifndef GOBLIN_FILM_H
#define GOBLIN_FILM_H

#include "GoblinColor.h"
#include <string>
namespace Goblin {
    class CameraSample;
    class Film {
    public:
        Film(int xRes, int yRes, const float crop[4], 
            const std::string& filename);
        ~Film(); 

        int getXResolution() const;
        int getYResolution() const;
        int getXStart() const;
        int getYStart() const;
        int getXEnd() const;
        int getYEnd() const;

        void addSample(const CameraSample& sample, const Color& L);
        void writeImage();

    private:
        int mXRes, mYRes;
        int mXStart, mYStart, mXEnd, mYEnd, mXCount, mYCount;
        float mCrop[4];
        std::string mFilename;
        class Pixel {
        public:
            Pixel(): color(0.0f, 0.0f, 0.0f, 1.0f), weight(0.0f) {}
            Color color;
            float weight;
            float pad[3];
        };
        Pixel* mPixels;
    };

    inline int Film::getXResolution() const { return mXRes; }
    inline int Film::getYResolution() const { return mYRes; }
    inline int Film::getXStart() const { return mXStart; }
    inline int Film::getXEnd() const { return mXEnd; }
    inline int Film::getYStart() const { return mYStart; }
    inline int Film::getYEnd() const { return mYEnd; }
}

#endif //GOBLIN_FILM_H
