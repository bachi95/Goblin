#ifndef GOBLIN_COLOR_H
#define GOBLIN_COLOR_H

#include "GoblinUtils.h"

#include <cassert>
#include <cmath>
#include <iostream>

namespace Goblin {

class Color {
public:
	static const Color Red;
	static const Color Green;
	static const Color Blue;
	static const Color White;
	static const Color Black;
	static const Color Yellow;
	static const Color Cyan;
	static const Color Magenta;
	float r, g, b, a;

public:
	Color() = default;

	explicit Color(float c) :r(c), g(c), b(c), a(1.0f) {}

	Color(float r, float g, float b, float a = 1.0f) :
		r(r), g(g), b(b), a(a) {}

	Color operator+(const Color& rhs) const {
		return Color(r + rhs.r, g + rhs.g, b + rhs.b, a);
	}

	Color& operator+=(const Color& rhs) {
		r += rhs.r;
		g += rhs.g;
		b += rhs.b;
		return *this;
	}

	Color operator-(const Color& rhs) const {
		return Color(r - rhs.r, g - rhs.g, b - rhs.b, a);
	}

	Color& operator-=(const Color& rhs) {
		r -= rhs.r;
		g -= rhs.g;
		b -= rhs.b;
		return *this;
	}

	Color operator*(float s) const {
		return Color(r * s, g * s, b * s, a);
	}

	Color& operator*=(float s) {
		r *= s;
		g *= s;
		b *= s;
		return *this;
	}

	Color operator*(const Color& rhs) const {
		return Color(r * rhs.r, g * rhs.g, b * rhs.b, a);
	}

	Color& operator*=(const Color& rhs) {
		r *= rhs.r;
		g *= rhs.g;
		b *= rhs.b;
		return *this;
	}

	Color operator/(float s) const {
		float inv = 1.0f / s;
		return Color(r * inv, g * inv, b * inv, a);
	}

	Color operator/(const Color& rhs) const {
		return Color(r / rhs.r, g / rhs.g, b / rhs.b, a);
	}

	Color& operator/=(const Color& rhs) {
		r /= rhs.r;
		g /= rhs.g;
		b /= rhs.b;
		return *this;
	}

	Color& operator/=(float s) {
		float inv = 1.0f / s;
		r *= inv;
		g *= inv;
		b *= inv;
		return *this;
	}

	Color operator-() const {
		return Color(-r, -g, -b, a);
	}

	bool operator==(const Color& rhs) const {
		return r == rhs.r && g == rhs.g && b == rhs.b && a == rhs.a;
	}

	bool operator!=(const Color& rhs) const {
		return !operator==(rhs);
	}

	bool isNaN() const {
		return Goblin::isNaN(r) ||
			Goblin::isNaN(g) ||
			Goblin::isNaN(b) ||
			Goblin::isNaN(a);
	}

	float luminance() const {
		return 0.212671f * r + 0.715160f * g + 0.072169f * b;
	}
};

inline Color operator*(float s, const Color& rhs) {
    return rhs * s;
}

inline Color sqrtColor(const Color& c) {
    return Color(sqrt(c.r), sqrt(c.g), sqrt(c.b));
}

inline Color expColor(const Color& c) {
    return Color(exp(c.r), exp(c.g), exp(c.b));
}

inline Color clampColor(const Color& c, float min = 0.0f,
    float max = INFINITY) {
    return Color(clamp(c.r, min, max),
        clamp(c.g, min, max),
        clamp(c.b, min, max));
}

inline std::ostream& operator<<(std::ostream& os, const Color& c) {
    os << "Color(" << c.r << ", " << c.g << ", " << c.b <<
        ", " << c.a << ")";
    return os;
}

}

#endif //GOBLIN_COLOR_H