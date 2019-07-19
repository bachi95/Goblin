#include "GoblinRenderContext.h"
#include "GoblinContextLoader.h"
#include <ctime>

using namespace Goblin;

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cout << "Usage: g_ray scene.json" << std::endl;
        return 0;
    }
    std::unique_ptr<RenderContext> renderContext(
        ContextLoader().load(argv[1]));
    if (renderContext) {
		std::cout << "\nsuccessfully loaded scene, start rendering..." <<
			std::endl;
        time_t beforeRender;
        time(&beforeRender);
        renderContext->render();
        time_t afterRender;
        time(&afterRender);
        double seconds = difftime(afterRender, beforeRender);
        std::cout << "render complete in " << seconds << " seconds!" <<
			std::endl;
    }
    return 0;
}
