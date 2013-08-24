#ifndef GOBLIN_IMAGEIO_H
#define GOBLIN_IMAGEIO_H

#include <string>

namespace Goblin {
    class Color;
    Color* loadImage(const std::string& filename, int *width, int *height);

    bool writeImage(const std::string& filename, const Color* colorBuffer,
        int width, int height);

}

#endif //GOBLIN_IMAGEIO_H 
