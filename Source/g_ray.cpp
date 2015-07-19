#include "GoblinRenderContext.h"
#include "GoblinContextLoader.h"

using namespace Goblin;

int main(int argc, char** argv) {
    if(argc != 2) {
        cout << "Usage: g_ray scene.json" << endl;
        return 0;
    }
    boost::scoped_ptr<RenderContext> renderContext(
        ContextLoader().load(argv[1]));
    if(renderContext) {
        cout << "\nsuccessfully loaded scene, start rendering...\n"; 
        renderContext->render();
        cout << "render complete!" << endl; 
    }
    return 0;
}
