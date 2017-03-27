#include "GoblinRenderContext.h"
#include "GoblinContextLoader.h"
#include <ctime>

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
        time_t beforeRender;
        time(&beforeRender);
        renderContext->render();
        time_t afterRender;
        time(&afterRender);
        double seconds = difftime(afterRender, beforeRender);
        cout << "render complete in " << seconds << " seconds!" << endl;
    }
    return 0;
}
