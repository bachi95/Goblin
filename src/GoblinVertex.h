#ifndef GOBLIN_VERTEX_H
#define GOBLIN_VERTEX_H

#include "GoblinVector.h"
namespace Goblin
{

struct Vertex {
    Vector3 position;
    Vector3 tangent;
    Vector3 normal;
    Vector2 texC;
};

}

#endif // GOBLIN_VERTEX_H