#ifndef GOBLIN_CONTEXT_LOADER_H
#define GOBLIN_CONTEXT_LOADER_H

#include <string>

namespace Goblin {
class RenderContext;

class ContextLoader {
public:
    static RenderContext* load(const std::string& filename);
};
}

#endif //GOBLIN_CONTEXT_LOADER_H
