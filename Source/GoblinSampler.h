#ifndef GOBLIN_SAMPLER_H
#define GOBLIN_SAMPLER_H

namespace Goblin {
    class CameraSample {
    public:
        CameraSample(float x, float y);
        // store film sample in image space( not NDC space)
        float imageX, imageY;
    };

    inline CameraSample::CameraSample(float x, float y):
        imageX(x), imageY(y) {}
}

#endif //GOBLIN_SAMPLER_H
