#include "GoblinGeometry.h"

namespace Goblin {
    size_t Geometry::nextGeometryId = 0;
    
    Geometry::Geometry(): mGeometryId(nextGeometryId++) {}

    const size_t Geometry::getId() const { return mGeometryId; }
}
