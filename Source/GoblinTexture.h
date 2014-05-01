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

    template<typename T>
    struct ImageBuffer {
        ImageBuffer(T* i, int w, int h): image(i), width(w), height(h) {}
        ~ImageBuffer() { delete [] image; image = NULL; }
        T texel(int s, int t, AddressMode addressMode);
        T lookup(float s, float t, AddressMode addressMode);
        T* image;
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

    template<typename T>
    class Texture {
    public:
        virtual ~Texture() {}
        virtual T lookup(const Fragment& f) const = 0;
    };

    typedef boost::shared_ptr<Texture<Color> > ColorTexturePtr;
    typedef boost::shared_ptr<Texture<float> > FloatTexturePtr;


    template<typename T>
    class ConstantTexture : public Texture<T> {
    public:
        ConstantTexture(const T& c);
        T lookup(const Fragment& f) const;
    private:
        T mValue;
    };

    template<typename T>
    class ScaleTexture : public Texture<T> {
    public:
        ScaleTexture(const boost::shared_ptr<Texture<T> >& t, 
            const FloatTexturePtr& s);
        T lookup(const Fragment& f) const;
    private:
        boost::shared_ptr<Texture<T> > mTexture;
        FloatTexturePtr mScale;
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

    template<typename T>
    class ImageTexture : public Texture<T> {
    public:
        ImageTexture(const string& filename, TextureMapping* m, 
            AddressMode address= AddressRepeat, float gamma = 1.0f);
        ~ImageTexture();
        T lookup(const Fragment& f) const;
        static void clearImageCache();
    private:
        ImageBuffer<T>* getImageBuffer(const TextureId& id);
    private:
        static std::map<TextureId, ImageBuffer<T>* > imageCache;
        TextureMapping* mMapping;
        AddressMode mAddressMode;
        ImageBuffer<T>* mImageBuffer;
    };

    template<typename T>
    void ImageTexture<T>::clearImageCache() {
        typename std::map<TextureId, ImageBuffer<T>* >::iterator it;
        for(it = imageCache.begin(); it != imageCache.end(); ++it) {
            std::cout << "clear image cache: " << it->first.filename << 
                std::endl;
            delete it->second;
        }
        imageCache.clear();
    }


    void gammaCorrect(const Color& in, float* out, float gamma);
    void gammaCorrect(const Color& in, Color* out, float gamma);

}

#endif //GOBLIN_TEXTURE_H
