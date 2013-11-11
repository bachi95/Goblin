#ifndef GOBLIN_PARAM_SET_H
#define GOBLIN_PARAM_SET_H
#include "GoblinColor.h"
#include "GoblinVector.h"

#include <string>
#include <vector>

/*
 * Emulation for python dict with limited types
 * for kwargs passing
 */

namespace Goblin {
    template <typename T> class ParamSetItem {
    public:
        ParamSetItem(const std::string& k, const T& d);
        std::string key;
        T data;
    };

    template <typename T>
    ParamSetItem<T>::ParamSetItem(const std::string& k, const T& d):
        key(k), data(d) {}

    class ParamSet {
    public:
        void setInt(const std::string& key, int i);
        void setFloat(const std::string& key, float f);
        void setVector3(const std::string& key, const Vector3& v);
        void setVector4(const std::string& key, const Vector4& v);
        void setColor(const std::string& key, const Color& c);
        void setString(const std::string& key, const std::string& s);

        bool eraseInt(const std::string& key);
        bool eraseFloat(const std::string& key);
        bool eraseVector3(const std::string& key);
        bool eraseVector4(const std::string& key);
        bool eraseColor(const std::string& key);
        bool eraseString(const std::string& key);

        int getInt(const std::string& key, int d = 0);
        float getFloat(const std::string& key, float d = 0.0f);
        Vector3 getVector3(const std::string& key, 
            const Vector3& d = Vector3::Zero);
        Vector4 getVector4(const std::string& key,
            const Vector4& d = Vector4::Zero);
        Color getColor(const std::string& key,
            const Color& d = Color::White);
        std::string getString(const std::string& key,
            const std::string& = "");

    private:
        std::vector<ParamSetItem<int> > mInts;
        std::vector<ParamSetItem<float> > mFloats;
        std::vector<ParamSetItem<Vector3> > mVec3s;
        std::vector<ParamSetItem<Vector4> > mVec4s;
        std::vector<ParamSetItem<Color> > mColors;
        std::vector<ParamSetItem<std::string> > mStrings;
    };
}
#endif //GOBLIN_PARAM_SET_H
