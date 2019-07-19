#ifndef GOBLIN_VERTEX_H
#define GOBLIN_VERTEX_H

#include "GoblinVector.h"
namespace Goblin
{
struct Vertex {
    Vertex() {}
    Vertex(float x, float y, float z, float nx, float ny, float nz,
        float u, float v): 
        position(x, y, z), normal(nx, ny, nz), texC(u, v) {}
        
    Vertex(float x, float y, float z, float tx, float ty, float tz,
        float nx, float ny, float nz, float u, float v):
        position(x, y, z), tangent(tx, ty, tz), normal(nx, ny, nz),
        texC(u, v) {}

    Vector3 position;
    Vector3 tangent;
    Vector3 normal;
    Vector2 texC;
};

struct VertexPNT {
    VertexPNT() {}
    VertexPNT(float x, float y, float z, float nx, float ny, float nz,
        float u, float v): 
        position(x, y, z), normal(nx, ny, nz), texC(u, v) {}

    Vector3 position;
    Vector3 normal;
    Vector2 texC;
};
}

#endif // GOBLIN_VERTEX_H
