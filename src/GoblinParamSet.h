#ifndef GOBLIN_PARAM_SET_H
#define GOBLIN_PARAM_SET_H
#include "GoblinColor.h"
#include "GoblinVector.h"

#include <string>
#include <vector>

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
    void setBool(const std::string& key, bool b);

    void setInt(const std::string& key, int i);

    void setFloat(const std::string& key, float f);

    void setVector2(const std::string& key, const Vector2& v);

    void setVector3(const std::string& key, const Vector3& v);

    void setVector4(const std::string& key, const Vector4& v);

    void setString(const std::string& key, const std::string& s);

    bool hasBool(const std::string& key) const;

    bool hasInt(const std::string& key) const;

    bool hasFloat(const std::string& key) const;

    bool hasVector2(const std::string& key) const;

    bool hasVector3(const std::string& key) const;

    bool hasVector4(const std::string& key) const;

    bool hasString(const std::string& key) const;

    bool getBool(const std::string& key, bool d = false) const;

    int getInt(const std::string& key, int d = 0) const;

    float getFloat(const std::string& key, float d = 0.0f) const;

    Vector2 getVector2(const std::string& key,
        const Vector2& d = Vector2::Zero) const;

    Vector3 getVector3(const std::string& key, 
        const Vector3& d = Vector3::Zero) const;

    Vector4 getVector4(const std::string& key,
        const Vector4& d = Vector4::Zero) const;

    std::string getString(const std::string& key,
        const std::string& = "") const;

private:
    std::vector<ParamSetItem<bool> > mBools;
    std::vector<ParamSetItem<int> > mInts;
    std::vector<ParamSetItem<float> > mFloats;
    std::vector<ParamSetItem<Vector2> > mVec2s;
    std::vector<ParamSetItem<Vector3> > mVec3s;
    std::vector<ParamSetItem<Vector4> > mVec4s;
    std::vector<ParamSetItem<std::string> > mStrings;
};

}

#endif //GOBLIN_PARAM_SET_H