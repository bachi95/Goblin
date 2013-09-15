#ifndef GOBLIN_SAMPLER_H
#define GOBLIN_SAMPLER_H

namespace Goblin {
    class CameraSample {
    public:
        CameraSample() {};
        CameraSample(float x, float y);
        // store film sample in image space (not NDC space)
        float imageX, imageY;
    };

    inline CameraSample::CameraSample(float x, float y):
        imageX(x), imageY(y) {}

    class Sampler {
    public:
        Sampler(int xStart, int xEnd, int yStart, int yEnd, 
            int xPerPixel, int yPerPixel);
        int maxSamplesPerRequest() const;
        int requestSamples(CameraSample* samples);
    private:
        void stratifiedSample2D(CameraSample* samples);

    private:
        int mXStart, mYStart;
        int mXEnd, mYEnd;
        int mCurrentX, mCurrentY;
        int mXPerPixel, mYPerPixel;
        int mSamplesPerPixel;

    };
}

#endif //GOBLIN_SAMPLER_H
