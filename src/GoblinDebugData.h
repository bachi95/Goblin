#ifndef GOBLIN_DEBUG_DATA_H
#define GOBLIN_DEBUG_DATA_H

#include "GoblinColor.h"
#include "GoblinRay.h"
#include "GoblinVector.h"

namespace Goblin {

typedef std::pair<Vector2, Vector2> DebugLine;

class DebugData {
public:
    void addRay(const Ray& ray, Color c = Color::White);
    void addPoint(const Vector3& point, Color c = Color::White);
    const std::vector<std::pair<Ray, Color> >& getRays() const;
    const std::vector<std::pair<Vector3, Color> >& getPoints() const;
private:
    std::vector<std::pair<Ray, Color> > mRays;
    std::vector<std::pair<Vector3, Color> > mPoints;
};

inline void DebugData::addRay(const Ray& line, Color c) {
    mRays.push_back(std::pair<Ray, Color>(line, c));
}

inline void DebugData::addPoint(const Vector3& point, Color c) {
    mPoints.push_back(std::pair<Vector3, Color>(point, c));
}

inline const std::vector<std::pair<Ray, Color> >& DebugData::getRays() const {
    return mRays;
}

inline const std::vector<std::pair<Vector3, Color> >& DebugData::getPoints() const {
    return mPoints;
}

}

#endif //GOBLIN_DEBUG_DATA_H