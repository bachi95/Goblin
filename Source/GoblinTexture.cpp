#include "GoblinTexture.h"
#include "GoblinGeometry.h"
#include "GoblinImageIO.h"
#include <cassert>

namespace Goblin{

    template<typename T>
    T ImageBuffer<T>::texel(int s, int t, AddressMode addressMode) {
        switch(addressMode) {
        case AddressClamp: {
            s = clamp(s, 0, width - 1);
            t = clamp(s, 0, height - 1);
            break;
        }
        case AddressBorder: {
            if(s < 0 || t < 0 || s >= width || t >= height) {
                return T(0.0f);
            }
            break;
        }
        case AddressRepeat: default: {
            s = s % width;
            t = t % height;
            if(s < 0) {
                s += width;
            }
            if(t < 0) {
                t += height;
            }
            break;
        }
        }
        return image[t * width + s];
    }

    template<typename T>
    T ImageBuffer<T>::lookup(float s, float t, AddressMode addressMode) {
        float sRes = s * width - 0.5f;
        float tRes = t * height - 0.5f;
        int s0 = floorInt(sRes);
        float ds = sRes - (float)s0;
        int t0 = floorInt(tRes);
        float dt = tRes - (float)t0;
        return (1.0f - ds) * (1.0f - dt) * texel(s0, t0, addressMode) +
               (       ds) * (1.0f - dt) * texel(s0 + 1, t0, addressMode) +
               (1.0f - ds) * (       dt) * texel(s0, t0 + 1, addressMode) +
               (       ds) * (       dt) * texel(s0 + 1, t0 + 1, addressMode);
    }

    UVMapping::UVMapping(const Vector2& scale, const Vector2& offset):
        mScale(scale), mOffset(offset) {}

    void UVMapping::map(const Fragment& f, TextureCoordinate* tc) const {
        const Vector2& uv = f.getUV();
        tc->st[0] = mScale[0] * uv[0] + mOffset[0];
        tc->st[1] = mScale[1] * uv[1] + mOffset[1];
    }

    template<typename T>
    ConstantTexture<T>::ConstantTexture(const T& c): mValue(c) {}

    template<typename T>
    T ConstantTexture<T>::lookup(const Fragment& f) const {
        return mValue;
    }

    template<typename T>
    ScaleTexture<T>::ScaleTexture(const boost::shared_ptr<Texture<T> >& t, 
        const FloatTexturePtr& s): mTexture(t), mScale(s) {}

    template<typename T>
    T ScaleTexture<T>::lookup(const Fragment& f) const {
        return mScale->lookup(f) * mTexture->lookup(f);
    }

    template<typename T>
    std::map<TextureId, ImageBuffer<T>* > ImageTexture<T>::imageCache;

    template<typename T>
    ImageTexture<T>::ImageTexture(const string& filename, TextureMapping* m,
        AddressMode address, float gamma): 
        mMapping(m), mAddressMode(address) {
        TextureId id(filename, gamma);
        mImageBuffer = getImageBuffer(id);
    }

    template<typename T>
    ImageTexture<T>::~ImageTexture() {
        if(mMapping) {
            delete mMapping;
        }
        mMapping = NULL;
    }

    template<typename T>
    T ImageTexture<T>::lookup(const Fragment& f) const {
        TextureCoordinate tc;
        mMapping->map(f, &tc);
        float s = tc.st[0];
        float t = tc.st[1];
        return mImageBuffer->lookup(s, t, mAddressMode);
    }

    void gammaCorrect(const Color& in, float* out, float gamma) {
        *out = pow(in.luminance(), gamma);
    }

    void gammaCorrect(const Color& in, Color* out, float gamma) {
        (*out).r = gamma == 1.0f ? in.r : pow(in.r, gamma);
        (*out).g = gamma == 1.0f ? in.g : pow(in.g, gamma);
        (*out).b = gamma == 1.0f ? in.b : pow(in.b, gamma);
        (*out).a = in.a;
    }

    template<typename T>
    ImageBuffer<T>* ImageTexture<T>::getImageBuffer(const TextureId& id) {
        if(imageCache.find(id) != imageCache.end()) {
            return imageCache[id];
        }
        int width, height;
        Color* colorBuffer = Goblin::loadImage(id.filename, &width, &height);
        T* texelBuffer;
        if(!colorBuffer) {
            std::cerr << "error loading image file " << 
                id.filename << std::endl;
            texelBuffer = new T[1];
            T texel;
            gammaCorrect(Color::Magenta, &texel, id.gamma);
            texelBuffer[0] = texel;
            width =  height = 1;
        } else {
            texelBuffer = new T[width * height];
            // gamma correct it back to linear
            for(int i = 0; i < width * height; ++i) {
                gammaCorrect(colorBuffer[i], &texelBuffer[i], id.gamma);
            }
            delete [] colorBuffer;
        }
        ImageBuffer<T>* ret = new ImageBuffer<T>(texelBuffer, width, height); 
        imageCache[id] = ret;
        return ret;
    }

    template struct ImageBuffer<float>;
    template struct ImageBuffer<Color>; 
    template class Texture<float>;
    template class Texture<Color>;
    template class ConstantTexture<float>;
    template class ConstantTexture<Color>;
    template class ScaleTexture<float>;
    template class ScaleTexture<Color>;
    template class ImageTexture<float>;
    template class ImageTexture<Color>;
}
