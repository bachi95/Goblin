#include "GoblinFilm.h"
#include "GoblinUtils.h"

#include <cstring>
#include <algorithm>

#include <iostream>

namespace Goblin {
    Film::Film(int xRes, int yRes, const float crop[4],
        const std::string& filename): mXRes(xRes), mYRes(yRes),
        mFilename(filename) {

        memcpy(mCrop, crop, 4 * sizeof(float));

        mXStart = ceilInt(mXRes * mCrop[0]);
        mXCount = std::max(1, ceilInt(mXRes * mCrop[1])- mXStart);
        mXEnd = mXStart + mXCount - 1;

        mYStart = ceilInt(mYRes * mCrop[2]); 
        mYCount = std::max(1, ceilInt(mYRes * mCrop[3]) - mYStart);
        mYEnd = mYStart + mYCount - 1;

        pixels = new Pixel[mXCount * mYCount];
    } 
    
    Film::~Film() {
        if(pixels != NULL) {
            delete[] pixels;
        }
    }

    void Film::addSample(int x, int y, const Color& L) {
        if(x < mXStart || x > mXEnd || y < mYStart || y < mYEnd) {
            std::cerr<< "(" << x << ", " << y <<")" <<
                "out of the film crop window" << std::endl;
            return;
        }
        //TODO atomic add when this come to multi thread
        //TODO filter table resolve for more delicate weight
        int index = y * mXRes +x;
        pixels[index].color += L;
        pixels[index].weight += 1.0f; 
    } 

    void Film::writeImage() {
    }

}
