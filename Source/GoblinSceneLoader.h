#ifndef GOBLIN_SCENE_LOADER_H
#define GOBLIN_SCENE_LOADER_H

#include "GoblinScene.h"
#include <string>

namespace Goblin {
    class SceneLoader {
    public:
        bool load(const std::string& filename, Scene* scene);

    };
}

#endif //GOBLIN_SCENE_LOADER_H