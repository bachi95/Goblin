#ifndef GOBLIN_COLOR_H
#define GOBLIN_COLOR_H

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
        Color();
        Color(float r, float g, float b, float a = 1.0f);
        Color operator+(const Color& rhs) const;
        Color& operator+=(const Color& rhs);
        Color operator-(const Color& rhs) const;
        Color& operator-=(const Color& rhs);
        Color operator*(float s) const;
        Color& operator*=(float s);
        Color operator*(const Color& rhs) const;
        Color& operator*=(const Color& rhs);
        Color operator/(float s) const;
        Color operator/(const Color& rhs) const;
        Color& operator/=(const Color& rhs);
        Color& operator/=(float s);
        Color operator-() const;
        bool operator==(const Color& rhs) const;
        bool operator!=(const Color& rhs) const;
        const float* ptr() const;
    };

    inline Color::Color() {}

    inline Color::Color(float r, float g, float b, float a) : 
        r(r), g(g), b(b), a(a) {}

    inline Color Color::operator+(const Color& rhs) const {
        return Color(r + rhs.r, g + rhs.g, b + rhs.b, a);
    }

    inline Color& Color::operator+=(const Color& rhs) {
        r += rhs.r; 
        g += rhs.g;
        b += rhs.b;
        return *this;
    }

    inline Color Color::operator-(const Color& rhs) const {
        return Color(r - rhs.r, g - rhs.g, b - rhs.b, a);
    }

    inline Color& Color::operator-=(const Color& rhs) {
        r -= rhs.r;
        g -= rhs.g;
        b -= rhs.b;
        return *this;
    }

    inline Color Color::operator*(float s) const {
        return Color(r * s, g * s, b * s, a);
    }

    inline Color& Color::operator*=(float s) {
        r *= s;
        g *= s;
        b *= s;
        return *this;
    }

    inline Color Color::operator*(const Color& rhs) const {
        return Color(r * rhs.r, g * rhs.g , b * rhs.b, a);
    }

    inline Color& Color::operator*=(const Color& rhs) {
        r *= rhs.r;
        g *= rhs.g;
        b *= rhs.b;
        return *this;
    }

    inline Color Color::operator/(float s) const {
        float inv = 1.0f / s;
        return Color(r * inv, g * inv, b * inv, a);
    }

    inline Color& Color::operator/=(float s) {
        float inv = 1.0f / s;
        r *= inv;
        g *= inv;
        b *= inv;
        return *this;
    }

    inline Color Color::operator/(const Color& rhs) const {
        return Color(r / rhs.r, g / rhs.g , b / rhs.b, a);
    }

    inline Color& Color::operator/=(const Color& rhs) {
        r /= rhs.r;
        g /= rhs.g;
        b /= rhs.b;
        return *this;
    }

    inline Color Color::operator-() const {
        return Color(-r, -g, -b, a);
    }

    inline bool Color::operator==(const Color& rhs) const {
        return r == rhs.r && g == rhs.g && b == rhs.b && a == rhs.a;
    }

    inline bool Color::operator!=(const Color& rhs) const {
        return !operator==(rhs);
    }

    inline const float* Color::ptr() const {
        return &r;
    }

    inline Color operator*(float s, const Color& rhs) {
        return rhs * s;
    }

    inline std::ostream& operator<<(std::ostream& os, const Color& c) {
        os << "Color(" << c.r << ", " << c.g << ", " << c.b << 
            ", " << c.a << ")";
        return os;
    }
}
#endif //GOBLIN_COLOR_H
