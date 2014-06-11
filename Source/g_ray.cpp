#include "GoblinScene.h"
#include "GoblinSceneLoader.h"
#include "GoblinWhitted.h"
#include "GoblinPathtracer.h"

#include <iostream>
#include <string>

using namespace std;
using namespace Goblin;

int main(int argc, char** argv) {
    if(argc != 2) {
        std::cout << "Usage: Ray scene.json" << std::endl;
        return 0;
    }
    RenderSetting setting;
    ScenePtr scene;
    if(scene = SceneLoader().load(argv[1], &setting)) {
        std::cout << "\nsuccessfully loaded scene, start rendering...\n"; 
        PathTracer renderer(setting);
        renderer.render(scene);
        std::cout << "render complete!" << std::endl; 
    }
    return 0;
}
