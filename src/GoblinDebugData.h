#ifndef GOBLIN_DEBUG_DATA_H
#define GOBLIN_DEBUG_DATA_H

#include "GoblinColor.h"
#include "GoblinRay.h"
#include "GoblinVector.h"

namespace Goblin {

typedef std::pair<Vector2, Vector2> DebugLine;

class DebugData {
public:
	void addRay(const Ray& ray, Color c = Color::White) {
		mRays.push_back(std::pair<Ray, Color>(ray, c));
	}

	void addPoint(const Vector3& point, Color c = Color::White) {
		mPoints.push_back(std::pair<Vector3, Color>(point, c));
	}

	const std::vector<std::pair<Ray, Color> >& getRays() const {
		return mRays;
	}

	const std::vector<std::pair<Vector3, Color> >& getPoints() const {
		return mPoints;
	}

private:
    std::vector<std::pair<Ray, Color> > mRays;
    std::vector<std::pair<Vector3, Color> > mPoints;
};

}

#endif //GOBLIN_DEBUG_DATA_H