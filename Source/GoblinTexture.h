#ifndef GOBLIN_TEXTURE_H
#define GOBLIN_TEXTURE_H

#include "GoblinColor.h"
#include "GoblinUtils.h"
#include "GoblinVector.h"
#include "GoblinTransform.h"
#include "GoblinFactory.h"

namespace Goblin {
    class Fragment;

    enum AddressMode {
        AddressRepeat,
        AddressClamp,
        AddressBorder
    };

    enum ImageChannel {
        ChannelR,
        ChannelG,
        ChannelB,
        ChannelA,
        ChannelAll
    };

    enum FilterType {
        FilterNone,
        FilterBilinear,
        FilterTrilinear,
        FilterEWA
    };

    struct TextureCoordinate {
        Vector2 st;
        float dsdx, dtdx, dsdy, dtdy;
    };

    template<typename T>
    struct ImageBuffer {
        ImageBuffer(T* i, int w, int h): image(i), width(w), height(h) {}
        ~ImageBuffer() { delete [] image; image = NULL; }
        T texel(int s, int t, AddressMode addressMode) const;
        T* image;
        int width, height;
    };

    template<typename T>
    class MIPMap {
    public:
        MIPMap(T* image, int w, int h, float maxAniso = 10.0f);
        ~MIPMap();
        T lookup(const TextureCoordinate& tc, 
            FilterType f, AddressMode m) const;
        // leave this explicitely select level lookup here for now...
        // sinc IBL try to figure out the level in its own math instead
        // of by ray differential
        T lookup(int level, float s, float t, 
            AddressMode m = AddressRepeat) const;
        const ImageBuffer<T>* getImageBuffer(int level) const;
        int getLevelsNum() const;
        int getWidth() const;
        int getHeight() const;

    private:
        T lookupNearest(float s, float t, AddressMode m) const;
        T lookupBilinear(float s, float t, float width, AddressMode m) const;
        T lookupTrilinear(float s, float t, float width, AddressMode m) const;
        T lookupEWA(const TextureCoordinate& tc, AddressMode m) const;
        T EWA(int level, float s, float t, float A, float B, float C, 
            AddressMode m) const;
        static void initEWALut();
    private:
        int mLevelsNum;
        int mWidth, mHeight;
        float mMaxAnisotropy;
        vector<ImageBuffer<T>* > mPyramid;

        static const size_t EWA_LUT_SIZE = 128;
        static vector<float> EWALut;
    };

    template<typename T>
    const ImageBuffer<T>* MIPMap<T>::getImageBuffer(int level) const {
        level = clamp(level, 0, mLevelsNum - 1);
        return mPyramid[level];
    }

    template<typename T>
    int MIPMap<T>::getLevelsNum() const { return mLevelsNum; }

    template<typename T>
    int MIPMap<T>::getWidth() const { return mWidth; }

    template<typename T>
    int MIPMap<T>::getHeight() const { return mHeight; }

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


    class SphericalMapping : public TextureMapping {
    public:
        SphericalMapping(const Transform& toTex);
        void map(const Fragment& f, TextureCoordinate* tc) const;
    private:
        void pointToST(const Vector3& p, float* s, float* t) const;
    private:
        Transform mToTex;
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
    class CheckboardTexture : public Texture<T> {
    public:
        CheckboardTexture(TextureMapping* m,
            const boost::shared_ptr<Texture<T> >& T1,
            const boost::shared_ptr<Texture<T> >& T2,
            bool mFilter);
        ~CheckboardTexture();
        T lookup(const Fragment& f) const;
    private:
        TextureMapping* mMapping;
        boost::shared_ptr<Texture<T> > mT1, mT2;
        bool mFilter;
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
        TextureId(const string& f, float g, ImageChannel c, float maxAniso): 
            filename(f), gamma(g), channel(c), maxAnisotropy(maxAniso) {}
        string filename;
        float gamma;
        ImageChannel channel;
        float maxAnisotropy;
        bool operator<(const TextureId &rhs) const;
    };

    inline bool TextureId::operator<(const TextureId & rhs) const {
        if(gamma != rhs.gamma) {
            return gamma < rhs.gamma;
        } else if(channel != rhs.channel) {
            return channel < rhs.channel;
        } else if(maxAnisotropy != rhs.maxAnisotropy) {
            return maxAnisotropy < rhs.maxAnisotropy;
        }
        return filename < rhs.filename;
    }

    template<typename T>
    class ImageTexture : public Texture<T> {
    public:
        ImageTexture(const string& filename, TextureMapping* m, 
            FilterType filter = FilterNone,
            AddressMode address= AddressRepeat, float gamma = 1.0f,
            ImageChannel channel = ChannelAll, float maxAnisotropy = 10.0f);
        ~ImageTexture();
        T lookup(const Fragment& f) const;
        static void clearImageCache();
    private:
        MIPMap<T>* getMIPMap(const TextureId& id);
        void convertTexel(const Color& in, T* out, float gamma, 
            ImageChannel channel);

    private:
        static std::map<TextureId, MIPMap<T>* > imageCache;
        TextureMapping* mMapping;
        FilterType mFilter;
        AddressMode mAddressMode;
        MIPMap<T>* mMIPMap;
    };

    template<typename T>
    void ImageTexture<T>::clearImageCache() {
        typename std::map<TextureId, MIPMap<T>* >::iterator it;
        for(it = imageCache.begin(); it != imageCache.end(); ++it) {
            delete it->second;
        }
        imageCache.clear();
    }

    class ParamSet;
    class SceneCache;

    class FloatConstantTextureCreator : public 
        Creator<Texture<float> , const ParamSet&, const SceneCache&> {
    public:
        Texture<float>* create(const ParamSet& params, 
            const SceneCache& sceneCache) const;
    };


    class FloatCheckboardTextureCreator : public 
        Creator<Texture<float> , const ParamSet&, const SceneCache&> {
    public:
        Texture<float>* create(const ParamSet& params, 
            const SceneCache& sceneCache) const;
    };


    class FloatScaleTextureCreator : public 
        Creator<Texture<float> , const ParamSet&, const SceneCache&> {
    public:
        Texture<float>* create(const ParamSet& params, 
            const SceneCache& sceneCache) const;
    };


    class FloatImageTextureCreator : public 
        Creator<Texture<float> , const ParamSet&, const SceneCache&> {
    public:
        Texture<float>* create(const ParamSet& params, 
            const SceneCache&) const;
    };


    class ColorConstantTextureCreator : public 
        Creator<Texture<Color> , const ParamSet&, const SceneCache&> {
    public:
        Texture<Color>* create(const ParamSet& params, 
            const SceneCache& sceneCache) const;
    };


    class ColorCheckboardTextureCreator : public 
        Creator<Texture<Color> , const ParamSet&, const SceneCache&> {
    public:
        Texture<Color>* create(const ParamSet& params, 
            const SceneCache& sceneCache) const;
    };


    class ColorScaleTextureCreator : public 
        Creator<Texture<Color> , const ParamSet&, const SceneCache&> {
    public:
        Texture<Color>* create(const ParamSet& params, 
            const SceneCache& sceneCache) const;
    };


    class ColorImageTextureCreator : public 
        Creator<Texture<Color> , const ParamSet&, const SceneCache&> {
    public:
        Texture<Color>* create(const ParamSet& params, 
            const SceneCache&) const;
    };


    template<typename T>
    T* resizeImage(const T* srcBuffer, int srcWidth, int srcHeight, 
        int dstWidth, int dstHeight);

    template<typename T>
    void convertTexel(const Color& in, T* out, float gamma, 
        ImageChannel channel);

}

#endif //GOBLIN_TEXTURE_H
