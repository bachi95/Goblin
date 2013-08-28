#ifndef GOBLIN_GEOMETRY_H
#define GOBLIN_GEOMETRY_H

#include <cstdio>

namespace Goblin {
    class Geometry {
    public:
        Geometry();
        virtual ~Geometry() {};
        virtual bool init() = 0;
        virtual const size_t getVertexNum() const = 0;
        virtual const size_t getFaceNum() const = 0;
        virtual const void* getVertexPtr() const = 0;
        virtual const void* getFacePtr() const = 0;
        const size_t getId() const;
    protected:
        static size_t nextGeometryId;
        size_t mGeometryId;
    };
}

#endif //GOBLIN_GEOMETRY_H
