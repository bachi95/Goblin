#ifndef GOBLIN_SCENE_LOADER_H
#define GOBLIN_SCENE_LOADER_H

#include "GoblinScene.h"
#include <string>

namespace Goblin {
    struct RenderSetting;
    class SceneLoader {
    public:
        ScenePtr load(const std::string& filename, 
            RenderSetting* setting = NULL);
    };
}

#endif //GOBLIN_SCENE_LOADER_H