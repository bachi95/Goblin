#ifndef GOBLIN_TEXTURE_H
#define GOBLIN_TEXTURE_H

#include "GoblinColor.h"
#include "GoblinUtils.h"
#include "GoblinVector.h"

#include <map>

namespace Goblin {
    class Fragment;

    enum AddressMode {
        AddressRepeat,
        AddressClamp,
        AddressBorder
    };

    struct TextureCoordinate {
        Vector2 st;
    };

    struct ImageBuffer {
        ImageBuffer(Color* i, int w, int h): image(i), width(w), height(h) {}
        ~ImageBuffer() { delete [] image; image = NULL; }
        Color texel(int s, int t, AddressMode addressMode);
        Color lookup(float s, float t, AddressMode addressMode);
        Color* image;
        int width, height;
    };

    class TextureMapping {
    public:
        virtual ~TextureMapping() {}
        virtual void map(const Fragment& f, TextureCoordinate* tc) const = 0;
    };

    class UVMapping : public TextureMapping {
    public:
        UVMapping(const Vector2& scale, const Vector2& offset);
        void map(const Fragment& f, TextureCoordinate* tc) const;
    private:
        Vector2 mScale;
        Vector2 mOffset;
    };

    class Texture {
    public:
        virtual ~Texture() {}
        virtual Color lookup(const Fragment& f) const = 0;
    };

    typedef boost::shared_ptr<Texture> TexturePtr;

    class ConstantTexture : public Texture {
    public:
        ConstantTexture(const Color& c);
        Color lookup(const Fragment& f) const;
    private:
        Color mValue;
    };

    struct TextureId {
        TextureId(const string& f, float g): filename(f), gamma(g) {}
        string filename;
        float gamma;
        bool operator<(const TextureId &rhs) const;
    };

    inline bool TextureId::operator<(const TextureId & rhs) const {
        if(gamma != rhs.gamma) {
            return gamma < rhs.gamma;
        }
        return filename < rhs.filename;
    }

    class ImageTexture : public Texture {
    public:
        ImageTexture(const string& filename, TextureMapping* m, 
            AddressMode address= AddressRepeat, float gamma = 1.0f);
        ~ImageTexture();
        Color lookup(const Fragment& f) const;
        static void clearImageCache();
    private:
        ImageBuffer* getImageBuffer(const TextureId& id);
    private:
        static std::map<TextureId, ImageBuffer*> imageCache;
        TextureMapping* mMapping;
        AddressMode mAddressMode;
        ImageBuffer* mImageBuffer;
    };
}

#endif //GOBLIN_TEXTURE_H
