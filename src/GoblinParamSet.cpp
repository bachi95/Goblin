#include "GoblinParamSet.h"

namespace Goblin {
void ParamSet::setBool(const std::string& key, bool b) {
    mBools.push_back(ParamSetItem<bool>(key, b));
}

void ParamSet::setInt(const std::string& key, int i) {
    mInts.push_back(ParamSetItem<int>(key, i));
}

void ParamSet::setFloat(const std::string& key, float f) {
    mFloats.push_back(ParamSetItem<float>(key, f));
}

void ParamSet::setVector2(const std::string& key, const Vector2& v) {
    mVec2s.push_back(ParamSetItem<Vector2>(key, v));
}

void ParamSet::setVector3(const std::string& key, const Vector3& v) {
    mVec3s.push_back(ParamSetItem<Vector3>(key, v));
}

void ParamSet::setVector4(const std::string& key, const Vector4& v) {
    mVec4s.push_back(ParamSetItem<Vector4>(key, v));
}

void ParamSet::setString(const std::string& key, const std::string& s) {
    mStrings.push_back(ParamSetItem<std::string>(key, s));
}

bool ParamSet::hasBool(const std::string& key) const {
    for (size_t i = 0; i < mBools.size(); ++i) {
        if (mBools[i].key == key) {
            return true;
        }
    }
    return false;
}

bool ParamSet::hasInt(const std::string& key) const {
    for (size_t i = 0; i < mInts.size(); ++i) {
        if (mInts[i].key == key) {
            return true;
        }
    }
    return false;
}

bool ParamSet::hasFloat(const std::string& key) const {
    for (size_t i = 0; i < mFloats.size(); ++i) {
        if (mFloats[i].key == key) {
            return true;
        }
    }
    return false;
}

bool ParamSet::hasVector2(const std::string& key) const {
    for (size_t i = 0; i < mVec2s.size(); ++i) {
        if (mVec2s[i].key == key) {
            return true;
        }
    }
    return false;
}

bool ParamSet::hasVector3(const std::string& key) const {
    for (size_t i = 0; i < mVec3s.size(); ++i) {
        if (mVec3s[i].key == key) {
            return true;
        }
    }
    return false;
}

bool ParamSet::hasVector4(const std::string& key) const {
    for (size_t i = 0; i < mVec4s.size(); ++i) {
        if (mVec4s[i].key == key) {
            return true;
        }
    }
    return false;
}

bool ParamSet::hasString(const std::string& key) const {
    for (size_t i = 0; i < mStrings.size(); ++i) {
        if (mStrings[i].key == key) {
            return true;
        }
    }
    return false;
}

bool ParamSet::getBool(const std::string& key, bool d) const {
    for (size_t i = 0; i < mBools.size(); ++i) {
        if (mBools[i].key == key) {
            return mBools[i].data;
        }
    }
    return d;
}

int ParamSet::getInt(const std::string& key, int d) const {
    for (size_t i = 0; i < mInts.size(); ++i) {
        if (mInts[i].key == key) {
            return mInts[i].data;
        }
    }
    return d;
}

float ParamSet::getFloat(const std::string& key, float d) const {
    for (size_t i = 0; i < mFloats.size(); ++i) {
        if (mFloats[i].key == key) {
            return mFloats[i].data;
        }
    }
    return d;
}

Vector2 ParamSet::getVector2(const std::string& key, 
    const Vector2& d) const {
    for (size_t i = 0; i < mVec2s.size(); ++i) {
        if (mVec2s[i].key == key) {
            return mVec2s[i].data;
        }
    }
    return d;
}

Vector3 ParamSet::getVector3(const std::string& key, 
    const Vector3& d) const {
    for (size_t i = 0; i < mVec3s.size(); ++i) {
        if (mVec3s[i].key == key) {
            return mVec3s[i].data;
        }
    }
    return d;
}

Vector4 ParamSet::getVector4(const std::string& key, 
    const Vector4& d) const {
    for (size_t i = 0; i < mVec4s.size(); ++i) {
        if (mVec4s[i].key == key) {
            return mVec4s[i].data;
        }
    }
    return d;
}

std::string ParamSet::getString(const std::string& key, 
    const std::string& d) const {
    for (size_t i = 0; i < mStrings.size(); ++i) {
        if (mStrings[i].key == key) {
            return mStrings[i].data;
        }
    }
    return d;
}
}
