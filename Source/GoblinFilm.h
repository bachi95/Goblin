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
}

#endif //GOBLIN_FILM_H
