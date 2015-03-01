#include "GoblinScene.h"
#include "GoblinParamSet.h"
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
        AddressMode address, float gamma, ImageChannel channel): 
        mMapping(m), mAddressMode(address) {
        TextureId id(filename, gamma, channel);
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
            convertTexel(Color::Magenta, &texel, id.gamma, id.channel);
            texelBuffer[0] = texel;
            width =  height = 1;
        } else {
            texelBuffer = new T[width * height];
            for(int i = 0; i < width * height; ++i) {
                convertTexel(colorBuffer[i], &texelBuffer[i], 
                    id.gamma, id.channel);
            }
            delete [] colorBuffer;
        }
        ImageBuffer<T>* ret = new ImageBuffer<T>(texelBuffer, width, height); 
        imageCache[id] = ret;
        return ret;
    }

    float gaussian(float x, float w, float falloff = 2.0f) {
        return max(0.0f, expf(-falloff * x * x) - expf(-falloff * w * w));
    }

    template<typename T>
    T* resizeImage(const T* srcBuffer, int srcWidth, int srcHeight,
        int dstWidth, int dstHeight) {
        float filterWidth = floor(max(2.0f, max(
            (float)srcWidth/(float)dstWidth,
            (float)srcHeight/(float)dstHeight)));
        int nSamples = floorInt(filterWidth) * 2;

        float* sSampleWeight = new float[dstWidth * nSamples];
        int* sSampleIndex = new int[dstWidth];
        for(int s = 0; s < dstWidth; ++s) {
            float center = ((float)s + 0.5f) / dstWidth * srcWidth;
            // the subtle from continuous space to discrete space offset...
            // if center at 2.75 (discrete space 2.25) you should sample
            // 1, 2, 3, 4, without that +0.5 offset would resolve to 0, 1, 2, 3
            sSampleIndex[s] = 
                floorInt(center - filterWidth + 0.5f);
            float weightSum = 0.0f;
            int wOffset = s * nSamples;
            for(int i = 0; i < nSamples; ++i) {
                float p = sSampleIndex[s] + 0.5f + i;
                sSampleWeight[wOffset + i] = gaussian(p - center, filterWidth);
                weightSum += sSampleWeight[wOffset + i];
            }
            float invW = 1.0f / weightSum;
            for(int i = 0; i < nSamples; ++i) {
                sSampleWeight[wOffset + i] *= invW;
            }
        }

        float* tSampleWeight = new float[dstHeight * nSamples];
        int* tSampleIndex = new int[dstHeight];
        for(int t = 0; t < dstHeight; ++t) {
            float center = ((float)t + 0.5f) / dstHeight * srcHeight;
            tSampleIndex[t] = 
                floorInt(center - filterWidth + 0.5f);
            float weightSum = 0.0f;
            int wOffset = t * nSamples;
            for(int i = 0; i < nSamples; ++i) {
                float p = tSampleIndex[t] + 0.5f + i;
                tSampleWeight[wOffset + i] = gaussian(p - center, filterWidth);
                weightSum += tSampleWeight[wOffset + i];
            }
            float invW = 1.0f / weightSum;
            for(int i = 0; i < nSamples; ++i) {
                tSampleWeight[wOffset + i] *= invW;
            }
        }

        T* dstBuffer = new T[dstWidth * dstHeight];
        for(int t = 0; t < dstHeight; ++t) {
            for(int s = 0; s < dstWidth; ++s) {
                int index = t * dstWidth + s;
                dstBuffer[index] = T(0.0f);
                for(int i = 0; i < nSamples; ++i) {
                    int srcT = clamp(tSampleIndex[t] + i, 0,
                        srcHeight - 1);
                    for(int j = 0; j < nSamples; ++j) {
                        int srcS = clamp(sSampleIndex[s] + j, 0,
                            srcWidth - 1);
                        float w = tSampleWeight[t * nSamples + i] * 
                            sSampleWeight[s * nSamples + j];
                        dstBuffer[index] += 
                            w * srcBuffer[srcT * srcWidth + srcS];
                    }
                }
            }
        }

        delete[] sSampleIndex;
        delete[] sSampleWeight;
        delete[] tSampleIndex;
        delete[] tSampleWeight;
        return dstBuffer;
    }

    static TextureMapping* getTextureMapping(const ParamSet& params) {
        TextureMapping* m;
        string type = params.getString("mapping");
        if(type == "uv") {
            Vector2 scale = params.getVector2("scale", Vector2(1.0f, 1.0f));
            Vector2 offset = params.getVector2("offset", Vector2::Zero);
            m = new UVMapping(scale, offset);
        } else {
            cerr << "undefined mapping type " << type << endl;
            m = new UVMapping(Vector2(1.0f, 1.0f), Vector2::Zero);
        }
        return m;
    }

    static AddressMode getAddressMode(const ParamSet& params) {
        string addressStr = params.getString("address", "repeat");
        AddressMode addressMode = AddressRepeat;
        if(addressStr == "repeat") {
            addressMode = AddressRepeat;
        } else if(addressStr == "clamp") {
            addressMode = AddressClamp;
        } else if(addressStr == "border") {
            addressMode = AddressBorder;
        }
        return addressMode;
    }

    Texture<float>* FloatConstantTextureCreator::create(const ParamSet& params,
            const SceneCache& sceneCache) const {
        float f = params.getFloat("float", 0.5f);
        return new ConstantTexture<float>(f);
    }

    Texture<float>* FloatScaleTextureCreator::create(const ParamSet& params,
            const SceneCache& sceneCache) const {
        string textureName = params.getString("texture");
        string scaleName = params.getString("scale");
        FloatTexturePtr s = sceneCache.getFloatTexture(scaleName);
        FloatTexturePtr t = sceneCache.getFloatTexture(textureName);
        return new ScaleTexture<float>(t, s);
    }

    static void getImageTextureParams(const ParamSet& params,
        const SceneCache& sceneCache,
        string* filePath, float* gamma, AddressMode* addressMode,
        ImageChannel* channel) {

        string filename = params.getString("file");
        *filePath = sceneCache.resolvePath(filename);
        *gamma = params.getFloat("gamma", 1.0f);
        string addressStr = params.getString("address", "repeat");
        *addressMode = getAddressMode(params);
        string channelStr = params.getString("channel", "All");
        *channel = ChannelAll;
        if(channelStr == "R") {
            *channel = ChannelR;
        } else if(channelStr == "G") {
            *channel = ChannelG;
        } else if(channelStr == "B") {
            *channel = ChannelB;
        } else if(channelStr == "A") {
            *channel = ChannelA;
        } else if(channelStr == "All") {
            *channel = ChannelAll;
        } else {
            cerr << "unrecognize channel: " << channelStr << endl;
        }
    }

    Texture<float>* FloatImageTextureCreator::create(const ParamSet& params,
            const SceneCache& sceneCache) const {
        TextureMapping* m = getTextureMapping(params);
        string filePath;
        float gamma;
        AddressMode addressMode;
        ImageChannel channel;
        getImageTextureParams(params, sceneCache, &filePath, &gamma, 
            &addressMode, &channel);
        return new ImageTexture<float>(filePath, m, 
            addressMode, gamma, channel);
    }

    Texture<Color>* ColorConstantTextureCreator::create(const ParamSet& params,
            const SceneCache& sceneCache) const {
        Color c = params.getColor("color", Color::Magenta);
        return new ConstantTexture<Color>(c);
    }

    Texture<Color>* ColorScaleTextureCreator::create(const ParamSet& params,
            const SceneCache& sceneCache) const {
        string textureName = params.getString("texture");
        string scaleName = params.getString("scale");
        FloatTexturePtr s = sceneCache.getFloatTexture(scaleName);
        ColorTexturePtr t = sceneCache.getColorTexture(textureName);
        return new ScaleTexture<Color>(t, s);
    }

    Texture<Color>* ColorImageTextureCreator::create(const ParamSet& params,
            const SceneCache& sceneCache) const {
        TextureMapping* m = getTextureMapping(params);
        string filePath;
        float gamma;
        AddressMode addressMode;
        ImageChannel channel;
        getImageTextureParams(params, sceneCache, &filePath, &gamma, 
            &addressMode, &channel);
        return new ImageTexture<Color>(filePath, m, 
            addressMode, gamma, channel);
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

    template Color* resizeImage<Color>(const Color* srcBuffer, 
        int srcWidth, int srcHeight,
        int dstWidth, int dstHeight);
    template float* resizeImage<float>(const float* srcBuffer, 
        int srcWidth, int srcHeight,
        int dstWidth, int dstHeight);

    template<>
    void convertTexel<float>(const Color& in, float* out, float gamma,
        ImageChannel channel) {
        if(channel == ChannelR) {
            *out = pow(in.r, gamma);
        } else if(channel == ChannelG) {
            *out = pow(in.g, gamma);
        } else if(channel == ChannelB) {
            *out = pow(in.b, gamma);
        } else if(channel == ChannelA) {
            *out = pow(in.a, gamma);
        } else if(channel == ChannelAll) {
            *out = pow(in.luminance(), gamma);
        }
    }

    template<>
    void convertTexel<Color>(const Color& in, Color* out, float gamma,
        ImageChannel channel) {
        Color c(in);
        if(channel == ChannelR) {
            c = Color(in.r);
        } else if(channel == ChannelG) {
            c = Color(in.g);
        } else if(channel == ChannelB) {
            c = Color(in.b);
        } else if(channel == ChannelA) {
            c = Color(in.a);
        }
        (*out).r = gamma == 1.0f ? c.r : pow(c.r, gamma);
        (*out).g = gamma == 1.0f ? c.g : pow(c.g, gamma);
        (*out).b = gamma == 1.0f ? c.b : pow(c.b, gamma);
        (*out).a = c.a;
    }
}
