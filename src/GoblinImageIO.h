#ifndef GOBLIN_IMAGEIO_H
#define GOBLIN_IMAGEIO_H

#include <string>

namespace Goblin {
class Color;
Color* loadImage(const std::string& filename, int *width, int *height);

bool writeImage(const std::string& filename, Color* colorBuffer,
    int width, int height, bool doToneMapping = false);

void bloom(Color* colorBuffer, int width, int height,
    float bloomRadius = 0.02f, float bloomWeight = 0.1f);

void toneMapping(Color* colorBuffer, int width, int height);
}

#endif //GOBLIN_IMAGEIO_H 