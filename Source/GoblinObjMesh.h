#ifndef GOBLIN_OBJ_MESH
#define GOBLIN_OBJ_MESH

#include "GoblinGeometry.h"
#include <string>

namespace Goblin {

    class ObjMesh : public Geometry {
    public:
        ObjMesh(const std::string& filename);
        ~ObjMesh();
        void init();
        bool load();

    private:
        std::string mFilename;
    };
}

#endif //GOBLIN_OBJ_MESH