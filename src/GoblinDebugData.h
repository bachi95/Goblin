#ifndef GOBLIN_DEBUG_DATA_H
#define GOBLIN_DEBUG_DATA_H

#include "GoblinColor.h"
#include "GoblinRay.h"
#include "GoblinVector.h"

namespace Goblin {

typedef pair<Vector2, Vector2> DebugLine;

class DebugData {
public:
    void addRay(const Ray& ray, Color c = Color::White);
    void addPoint(const Vector3& point, Color c = Color::White);
    const vector<pair<Ray, Color> >& getRays() const;
    const vector<pair<Vector3, Color> >& getPoints() const;
private:
    vector<pair<Ray, Color> > mRays;
    vector<pair<Vector3, Color> > mPoints;
};

inline void DebugData::addRay(const Ray& line, Color c) {
    mRays.push_back(pair<Ray, Color>(line, c));
}

inline void DebugData::addPoint(const Vector3& point, Color c) {
    mPoints.push_back(pair<Vector3, Color>(point, c));
}

inline const vector<pair<Ray, Color> >& DebugData::getRays() const {
    return mRays;
}

inline const vector<pair<Vector3, Color> >& DebugData::getPoints() const {
    return mPoints;
}

}

#endif //GOBLIN_DEBUG_DATA_H