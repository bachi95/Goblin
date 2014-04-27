#include "GoblinTexture.h"
#include "GoblinGeometry.h"
#include "GoblinImageIO.h"
#include <cassert>

namespace Goblin{
    Color ImageBuffer::texel(int s, int t, AddressMode addressMode) {
        switch(addressMode) {
        case AddressClamp: {
            s = clamp(s, 0, width - 1);
            t = clamp(s, 0, height - 1);
            break;
        }
        case AddressBorder: {
            if(s < 0 || t < 0 || s >= width || t >= height) {
                return Color::Black;
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

    Color ImageBuffer::lookup(float s, float t, AddressMode addressMode) {
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

    ConstantTexture::ConstantTexture(const Color& c): mValue(c) {}

    Color ConstantTexture::lookup(const Fragment& f) const {
        return mValue;
    }

    std::map<TextureId, ImageBuffer*> ImageTexture::imageCache;

    ImageTexture::ImageTexture(const string& filename, TextureMapping* m,
        AddressMode address, float gamma): 
        mMapping(m), mAddressMode(address) {
        TextureId id(filename, gamma);
        mImageBuffer = getImageBuffer(id);
    }

    ImageTexture::~ImageTexture() {
        if(mMapping) {
            delete mMapping;
        }
        mMapping = NULL;
    }

    Color ImageTexture::lookup(const Fragment& f) const {
        TextureCoordinate tc;
        mMapping->map(f, &tc);
        float s = tc.st[0];
        float t = tc.st[1];
        return mImageBuffer->lookup(s, t, mAddressMode);
    }

    void ImageTexture::clearImageCache() {
        std::map<TextureId, ImageBuffer*>::iterator it;
        for(it = imageCache.begin(); it != imageCache.end(); ++it) {
            std::cout << "clear image cache: " << it->first.filename << 
                std::endl;
            delete it->second;
        }
        imageCache.clear();
    }

    ImageBuffer* ImageTexture::getImageBuffer(const TextureId& id) {
        if(imageCache.find(id) != imageCache.end()) {
            return imageCache[id];
        }
        int width, height;
        Color* colorBuffer = Goblin::loadImage(id.filename, &width, &height);
        if(!colorBuffer) {
            std::cerr << "error loading image file " << 
                id.filename << std::endl;
            colorBuffer = new Color[1];
            colorBuffer[0] = Color::Magenta;
            width =  height = 1;
        }

        if(id.gamma != 1.0f) {
            // gamma correct it back to linear
            for(int i = 0; i < width * height; ++i) {
                colorBuffer[i].r = pow(colorBuffer[i].r, id.gamma);
                colorBuffer[i].g = pow(colorBuffer[i].g, id.gamma);
                colorBuffer[i].b = pow(colorBuffer[i].b, id.gamma);
            }
        }
        ImageBuffer* ret = new ImageBuffer(colorBuffer, width, height); 
        imageCache[id] = ret;
        return ret;
    }
}