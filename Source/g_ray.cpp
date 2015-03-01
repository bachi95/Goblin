#include "GoblinScene.h"
#include "GoblinSceneLoader.h"
#include "GoblinWhitted.h"
#include "GoblinPathtracer.h"

using namespace Goblin;

int main(int argc, char** argv) {
    if(argc != 2) {
        cout << "Usage: Ray scene.json" << endl;
        return 0;
    }
    ParamSet setting;
    ScenePtr scene;
    if(scene = SceneLoader().load(argv[1], &setting)) {
        cout << "\nsuccessfully loaded scene, start rendering...\n"; 
        // TODO make this a factory method when we have more renderers...
        string method = setting.getString("render_method", "path_tracing");
        Renderer* renderer;
        if(method == "whitted") {
            renderer = new WhittedRenderer(setting);
        } else {
            renderer = new PathTracer(setting);
        }
        renderer->render(scene);
        cout << "render complete!" << endl; 
        delete renderer;
    }
    return 0;
}
