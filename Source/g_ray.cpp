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
        cout << "Usage: Ray scene.json" << endl;
        return 0;
    }
    RenderSetting setting;
    ScenePtr scene;
    if(scene = SceneLoader().load(argv[1], &setting)) {
        cout << "\nsuccessfully loaded scene, start rendering...\n"; 
        // TODO make this a factory method when we have more renderers...
        Renderer* renderer;
        if(setting.method == PathTracing) {
            renderer = new PathTracer(setting);
        } else {
            renderer = new WhittedRenderer(setting);
        }
        renderer->render(scene);
        cout << "render complete!" << endl; 
        delete renderer;
    }
    return 0;
}
