#include "GoblinScene.h"
#include "GoblinParamSet.h"
#include "GoblinTexture.h"
#include "GoblinGeometry.h"
#include "GoblinImageIO.h"
#include <cassert>

namespace Goblin{

    template<typename T>
    T ImageBuffer<T>::texel(int s, int t, AddressMode addressMode) const {
        switch(addressMode) {
        case AddressClamp: {
            s = clamp(s, 0, width - 1);
            t = clamp(s, 0, height - 1);
            break;
        }
        case AddressBorder: {
            if (s < 0 || t < 0 || s >= width || t >= height) {
                return T(0.0f);
            }
            break;
        }
        case AddressRepeat: default: {
            s = s % width;
            t = t % height;
            if (s < 0) {
                s += width;
            }
            if (t < 0) {
                t += height;
            }
            break;
        }
        }
        return image[t * width + s];
    }

    template<typename T>
    MIPMap<T>::MIPMap(T* image, int w, int h, float maxAniso): 
        mWidth(w), mHeight(h), mMaxAnisotropy(maxAniso) {
        // resize the image to power of 2 for easier mipmap process
        if (!isPowerOf2(mWidth) || !isPowerOf2(mHeight)) {
            int wPow2 = roundUpPow2(mWidth);
            int hPow2 = roundUpPow2(mHeight);
            T* resizeBuffer = 
                resizeImage(image, mWidth, mHeight, wPow2, hPow2);
            delete [] image;
            image = resizeBuffer;
            mWidth = wPow2;
            mHeight = hPow2;
        }
        // start building the mipmap
        mPyramid.push_back(new ImageBuffer<T>(image, mWidth, mHeight));
        mLevelsNum = floorInt(
            max(log2((float)mWidth), log2((float)mHeight))) + 1;
        for (int i = 1; i < mLevelsNum; ++i) {
            const ImageBuffer<T>* buffer = mPyramid[i - 1];
            int resizedW = max(1, buffer->width >> 1);
            int resizedH = max(1, buffer->height >> 1);
            T* resizedBuffer = resizeImage(buffer->image, 
                buffer->width, buffer->height,
                resizedW, resizedH);
            mPyramid.push_back(
                new ImageBuffer<T>(resizedBuffer, resizedW, resizedH));
        }

        if (EWALut.empty()) {
            initEWALut();
        }
    }

    template<typename T>
    MIPMap<T>::~MIPMap() {
        for (size_t i = 0; i < mPyramid.size(); ++i) {
            delete mPyramid[i];
            mPyramid[i] = NULL;
        }
        mPyramid.clear();
    }

    template<typename T>
    T MIPMap<T>::lookup(const TextureCoordinate& tc, 
        FilterType f, AddressMode m) const {
        if (f == FilterNone) {
            return lookupNearest(tc.st[0], tc.st[1], m);
        } else if (f == FilterBilinear) {
            float width = max(max(fabs(tc.dsdx), fabs(tc.dtdx)), 
                max(fabs(tc.dsdy), fabs(tc.dtdy)));
            return lookupBilinear(tc.st[0], tc.st[1], width, m);
        } else if (f == FilterTrilinear) {
            float width = max(max(fabs(tc.dsdx), fabs(tc.dtdx)), 
                max(fabs(tc.dsdy), fabs(tc.dtdy)));
            return lookupTrilinear(tc.st[0], tc.st[1], width, m);
        } else if (f == FilterEWA) {
            return lookupEWA(tc, m);
        } else {
            return lookupNearest(tc.st[0], tc.st[1], m);
        }
    }

    template<typename T>
    T MIPMap<T>::lookupNearest(float s, float t, AddressMode m) const {
        return lookup(0, s, t, m);
    }

    template<typename T>
    T MIPMap<T>::lookupBilinear(float s, float t, float width,
        AddressMode m) const {
        float level = mLevelsNum - 1 + log2(max(width, 1e-8f));
        int iLevel = roundInt(level);
        return lookup(iLevel, s, t, m);
    }

    template<typename T>
    T MIPMap<T>::lookupTrilinear(float s, float t, float width,
        AddressMode m) const {
        float level = mLevelsNum - 1 + log2(max(width, 1e-8f));
        int iLevel = floorInt(level);
        if (iLevel < 0) {
            return lookup(0, s, t, m);
        } else if (iLevel >= mLevelsNum - 1) {
            return lookup(mLevelsNum - 1, s, t, m);
        } else {
            float delta = level - (float)iLevel;
            return (1.0f - delta) * lookup(iLevel, s, t, m) +
                   (       delta) * lookup(iLevel + 1, s, t, m);
        }
    }

    template<typename T>
    T MIPMap<T>::lookupEWA(const TextureCoordinate& tc, AddressMode m) const {
        float ds0 = tc.dsdx;
        float dt0 = tc.dtdx;
        float ds1 = tc.dsdy;
        float dt1 = tc.dtdy;
        float majorLength = sqrt(ds0 * ds0 + dt0 * dt0);
        float minorLength = sqrt(ds1 * ds1 + dt1 * dt1);
        if (majorLength < minorLength) {
            swap(ds0, ds1);
            swap(dt0, dt1);
            swap(majorLength, minorLength);
        }
        // clamp the eccentricity to a reasonable number
        if (minorLength * mMaxAnisotropy < majorLength && minorLength > 0.0f) {
            float scale = majorLength / (minorLength * mMaxAnisotropy);
            minorLength *= scale;
            ds1 *= scale;
            dt1 *= scale;
        }
        /* 
         * compute the ellipse equation coefficient
         * ellipse equation Ax^2 + Bxy +Cy^2 = F in matrix form:
         * P = [x, y]
         * Q = [A    B/2]
         *     [B/2  C  ]
         * PQPt = F (t stands for transpose)
         * now mapping a unit sphere to space formed by (ds0, dt0), (ds1, dt1)
         * ie: mapping (1, 0) to (ds0, dt0)
         *             (0, 1) to (ds1, dt1)
         * P' = PM = [cosTheta, sinTheta] [ds0 dt0]
         *                                [ds1 dt1]
         * where M is the transform matrix (jacobian)
         * since P = P'M^-1 and PPt = [cosTheta sinTheta] [cosTheta] = 1
         *                                                [sinTheta]
         * P'M^-1M^1tP't = 1 -> M^-1M^-1t = Q
         * Q = [A    B/2] = [ dt1 -dt0][ dt1 -ds1] ->
         *     [B/2  C  ]   [-ds1  ds0][-dt0  ds0] 
         * A = dt0^2 + dt1^2
         * B = -2(ds0dt0 + ds1dt1)
         * C = ds0^2 + ds1^2
         * F = (ds0dt1 -dt0ds1)^2 = AC - B^2/4
         */
        float A = dt0 * dt0 + dt1 * dt1;
        float B = -2.0f * (ds0 * dt0 + ds1 * dt1);
        float C = ds0 * ds0 + ds1 * ds1;
        float F =  A * C - 0.25f * B * B;
        // unable to form an ellipse, fallback to trilinear filter then
        float s = tc.st[0];
        float t = tc.st[1];
        if (minorLength == 0.0f || F <= 0.0f) {
            return lookupTrilinear(s, t, minorLength, m);
        }
        float invF = 1.0f / F;
        A *= invF;
        B *= invF;
        C *= invF;
        float level = mLevelsNum - 1 + log2(minorLength);
        int iLevel = floorInt(level);
        if (iLevel < 0) {
            return lookup(0, s, t, m);
        } else if (iLevel >= mLevelsNum - 1) {
            return lookup(mLevelsNum - 1, s, t, m);
        } else {
            float delta = level - (float)iLevel;
            return (1.0f - delta) * EWA(iLevel, s, t, A, B, C, m) +
                   (       delta) * EWA(iLevel + 1, s, t, A, B, C, m);
        }
    }

    template<typename T>
    T MIPMap<T>::EWA(int level, float s, float t, float A, float B, float C,
        AddressMode m) const {
        const ImageBuffer<T>* image = mPyramid[level];
        float sRes = (float)image->width;
        float tRes = (float)image->height;
        s = s * image->width - 0.5f;
        t = t * image->height - 0.5f;
        /*
         * transform the ellipse coefficient to pixel space
         * since dt0' = dt0 * tRes   ds0' = ds0 * sRes
         *       dt1' = dt1 * tRes   ds1' = ds1 * sRes
         * plug in back to the transform equation and we can get
         * (A'/F') = (A/F)/sRes^2
         * (B'/F') = (B/F)/(sRestRes)
         * (C'/F') = (C/F)/tRes^2
         */
        A = A / (sRes * sRes);
        B = B / (sRes * tRes);
        C = C / (tRes * tRes);
        /* 
         * compute the bounding box of ellipse
         * f = As^2 + Bst + Ct^2 -1
         * s bound happens in wehre dfdt = Bs + 2Ct = 0 ->
         * t = -Bs/2C -> As^2 + -B^2s^2/2C + B^2Cs^2/4C^2 - 1 = 0
         * s = +-2sqrt(C/(-B^2 + 4AC))
         * t bound happens in where dfds = 2As + Bt = 0 ->
         * s = -Bt/2A -> AB^2t^2/4A^2 + -B^2t^2/2A + Ct^2 - 1 = 0
         * t = +-2sqrt(A/(-B^2 + 4AC))
         */
        float invDet = 1.0f / (-B * B + 4.0f * A * C);
        float offsetS = 2.0f * sqrt(C * invDet);
        float offsetT = 2.0f * sqrt(A * invDet);
        int s0 = ceilInt(s - offsetS);
        int s1 = floorInt(s + offsetS);
        int t0 = ceilInt(t - offsetT);
        int t1 = floorInt(t + offsetT);
        // loop over bounding box region and accumulate pixel contribution
        float weightSum = 0.0f;
        T result(0.0f);
        for (int is = s0; is <= s1; ++is) {
            for (int it = t0; it <= t1; ++it) {
                float ss = is - s;
                float tt = it - t;
                float r2 = A * ss * ss + B * ss * tt + C * tt * tt;
                if (r2 <= 1.0f) {
                    size_t lutIndex = (size_t)floorInt(r2 * EWA_LUT_SIZE);
                    float weight = EWALut[min(lutIndex, EWA_LUT_SIZE - 1)];
                    result += weight * image->texel(is, it, m);
                    weightSum += weight;
                }
            }
        }
        if (weightSum > 0.0f) {
            result /= weightSum;
        } else {
            result = image->texel((int)s, (int)t, m);
        }
        return result;
    }

    template<typename T>
    void MIPMap<T>::initEWALut() {
        EWALut.resize(EWA_LUT_SIZE);
        const float falloff = 2.0f;
        for (size_t i = 0; i < EWA_LUT_SIZE; ++i) {
            float r2 = float(i) / float(EWA_LUT_SIZE - 1);
            // offset the whole lut so the r2 >= 1 cases have weight 0
            EWALut[i] = expf(-falloff * r2) - expf(-falloff);
        }
    }

    template<typename T>
    T MIPMap<T>::lookup(int level, float s, float t, 
        AddressMode m) const {
        level = clamp(level, 0, mLevelsNum);
        const ImageBuffer<T>* image = mPyramid[level];
        float sRes = s * image->width - 0.5f;
        float tRes = t * image->height - 0.5f;
        int s0 = floorInt(sRes);
        float ds = sRes - (float)s0;
        int t0 = floorInt(tRes);
        float dt = tRes - (float)t0;
        return (1.0f - ds) * (1.0f - dt) * image->texel(s0, t0, m) +
               (       ds) * (1.0f - dt) * image->texel(s0 + 1, t0, m) +
               (1.0f - ds) * (       dt) * image->texel(s0, t0 + 1, m) +
               (       ds) * (       dt) * image->texel(s0 + 1, t0 + 1, m);
    }

    template<typename T>
    std::vector<float> MIPMap<T>::EWALut;

    UVMapping::UVMapping(const Vector2& scale, const Vector2& offset):
        mScale(scale), mOffset(offset) {}

    void UVMapping::map(const Fragment& f, TextureCoordinate* tc) const {
        const Vector2& uv = f.getUV();
        tc->st[0] = mScale[0] * uv[0] + mOffset[0];
        tc->st[1] = mScale[1] * uv[1] + mOffset[1];
        tc->dsdx = mScale[0] * f.getDUDX();
        tc->dtdx = mScale[1] * f.getDVDX();
        tc->dsdy = mScale[0] * f.getDUDY();
        tc->dtdy = mScale[1] * f.getDVDY();
    }

    SphericalMapping::SphericalMapping(const Transform& toTex): 
        mToTex(toTex) {}

    void SphericalMapping::map(const Fragment& f, TextureCoordinate* tc) const {
        float s, t;
        const Vector3& p = f.getPosition();
        pointToST(p, &s, &t);
        tc->st[0] = s;
        tc->st[1] = t;
        // foward difference approximation, delta 1 pixel
        float sdx, tdx;
        pointToST(p + f.getDPDX(), &sdx, &tdx);
        float sdy, tdy;
        pointToST(p + f.getDPDY(), &sdy, &tdy);
        // take care of the discontinuous phi angle transition
        float dsdx = sdx - s;
        if (dsdx > 0.5f) {
            dsdx -= 1.0f;
        } else if (dsdx < -0.5f) {
            dsdx += 1.0f;
        }
        float dsdy = sdy -s;
        if (dsdy > 0.5f) {
            dsdy -= 1.0f;
        } else if (dsdy < -0.5f) {
            dsdy += 1.0f;
        }
        tc->dsdx = dsdx;
        tc->dtdx = tdx - t;
        tc->dsdy = dsdy;
        tc->dtdy = tdy - t;
    }

    void SphericalMapping::pointToST(const Vector3& p, 
        float* s, float* t) const {
        Vector3 v = normalize(mToTex.onPoint(p) - Vector3::Zero);
        float theta = sphericalTheta(v);
        float phi = sphericalPhi(v);
        *s = phi * INV_TWOPI;
        *t = theta * INV_PI;
    }

    template<typename T>
    ConstantTexture<T>::ConstantTexture(const T& c): mValue(c) {}

    template<typename T>
    T ConstantTexture<T>::lookup(const Fragment& f) const {
        return mValue;
    }

    template<typename T>
    CheckboardTexture<T>::CheckboardTexture(TextureMapping* m,
        const std::shared_ptr<Texture<T> >& T1,
        const std::shared_ptr<Texture<T> >& T2,
        bool filter):
        mMapping(m), mT1(T1), mT2(T2), mFilter(filter) {}

    template<typename T>
    CheckboardTexture<T>::~CheckboardTexture() {
        if (mMapping) {
            delete mMapping;
        }
        mMapping = NULL;
    }

    static float integrateChecker(float x) {
        float xHalf = 0.5f * x;
        return floor(xHalf) + 
            2.0f * max(xHalf - floor(xHalf) - 0.5f, 0.0f);
    }

    template<typename T>
    T CheckboardTexture<T>::lookup(const Fragment& f) const {
        TextureCoordinate tc;
        mMapping->map(f, &tc);
        float s = tc.st[0];
        float t = tc.st[1];
        if (!mFilter) {
            return (floorInt(s) + floorInt(t)) % 2 == 0 ?
                mT1->lookup(f) : mT2->lookup(f);
        }
        float ds = max(fabs(tc.dsdx), fabs(tc.dsdy));
        float dt = max(fabs(tc.dtdx), fabs(tc.dtdy));
        float s0 = s - ds;
        float s1 = s + ds;
        float t0 = t - dt;
        float t1 = t + dt;
        if (floorInt(s0) == floorInt(s1) && floorInt(t0) == floorInt(t1)) {
            return (floorInt(s) + floorInt(t)) % 2 == 0 ?
                mT1->lookup(f) : mT2->lookup(f);
        } else {
            float sTex2Ratio = (integrateChecker(s1) - integrateChecker(s0)) / 
                (2.0f * ds);
            float tTex2Ratio = (integrateChecker(t1) - integrateChecker(t0)) /
                (2.0f * dt);
            /*
             * draw a picture should make this clearer...but imagine this
             * ds * dt uv box fall on the intersection of checker boarder
             * the textur2 ratio of the total uv box would be:
             * sTex2Ratio * (1 - tTex2Ratio) + tTex2Ratio * (1 - sTex2Ratio) =
             * sTex2Ratio + tTex2Ratio - 2 * sTex2Ratio * tTex2Ratio
             */
            float tex2Area = sTex2Ratio + tTex2Ratio - 
                2.0f * sTex2Ratio * tTex2Ratio;
            if (ds > 1.0f || dt > 1.0f) {
                tex2Area = 0.5f;
            }
            return (1.0f - tex2Area) * mT1->lookup(f) + 
                tex2Area * mT2->lookup(f);
        }
    }

    template<typename T>
    ScaleTexture<T>::ScaleTexture(const std::shared_ptr<Texture<T> >& t, 
        const FloatTexturePtr& s): mTexture(t), mScale(s) {}

    template<typename T>
    T ScaleTexture<T>::lookup(const Fragment& f) const {
        return mScale->lookup(f) * mTexture->lookup(f);
    }

    template<typename T>
    std::map<TextureId, MIPMap<T>* > ImageTexture<T>::imageCache;

    template<typename T>
    ImageTexture<T>::ImageTexture(const string& filename, TextureMapping* m,
        FilterType filter, AddressMode address, 
        float gamma, ImageChannel channel, float maxAniso): 
        mMapping(m), mFilter(filter), mAddressMode(address) {
        TextureId id(filename, gamma, channel, maxAniso);
        mMIPMap = getMIPMap(id);
    }

    template<typename T>
    ImageTexture<T>::~ImageTexture() {
        if (mMapping) {
            delete mMapping;
        }
        mMapping = NULL;
    }

    template<typename T>
    T ImageTexture<T>::lookup(const Fragment& f) const {
        TextureCoordinate tc;
        mMapping->map(f, &tc);
        return mMIPMap->lookup(tc, mFilter, mAddressMode);
    }

    template<typename T>
    MIPMap<T>* ImageTexture<T>::getMIPMap(const TextureId& id) {
        if (imageCache.find(id) != imageCache.end()) {
            return imageCache[id];
        }
        int width, height;
        Color* colorBuffer = Goblin::loadImage(id.filename, &width, &height);
        T* texelBuffer;
        if (!colorBuffer) {
            std::cerr << "error loading image file " << 
                id.filename << std::endl;
            texelBuffer = new T[1];
            T texel;
            convertTexel(Color::Magenta, &texel, id.gamma, id.channel);
            texelBuffer[0] = texel;
            width =  height = 1;
        } else {
            texelBuffer = new T[width * height];
            for (int i = 0; i < width * height; ++i) {
                convertTexel(colorBuffer[i], &texelBuffer[i], 
                    id.gamma, id.channel);
            }
            delete [] colorBuffer;
        }
        MIPMap<T>* ret = new MIPMap<T>(texelBuffer, width, height, 
            id.maxAnisotropy);
        imageCache[id] = ret;
        return ret;
    }

    template<>
    void ImageTexture<float>::convertTexel(const Color& in, float* out, 
        float gamma, ImageChannel channel) {
        if (channel == ChannelR) {
            *out = pow(in.r, gamma);
        } else if (channel == ChannelG) {
            *out = pow(in.g, gamma);
        } else if (channel == ChannelB) {
            *out = pow(in.b, gamma);
        } else if (channel == ChannelA) {
            *out = pow(in.a, gamma);
        } else if (channel == ChannelAll) {
            *out = pow(in.luminance(), gamma);
        }
    }

    template<>
    void ImageTexture<Color>::convertTexel(const Color& in, Color* out, 
        float gamma, ImageChannel channel) {
        Color c(in);
        if (channel == ChannelR) {
            c = Color(in.r);
        } else if (channel == ChannelG) {
            c = Color(in.g);
        } else if (channel == ChannelB) {
            c = Color(in.b);
        } else if (channel == ChannelA) {
            c = Color(in.a);
        }
        (*out).r = gamma == 1.0f ? c.r : pow(c.r, gamma);
        (*out).g = gamma == 1.0f ? c.g : pow(c.g, gamma);
        (*out).b = gamma == 1.0f ? c.b : pow(c.b, gamma);
        (*out).a = c.a;
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
        for (int s = 0; s < dstWidth; ++s) {
            float center = ((float)s + 0.5f) / dstWidth * srcWidth;
            // the subtle from continuous space to discrete space offset...
            // if center at 2.75 (discrete space 2.25) you should sample
            // 1, 2, 3, 4, without that +0.5 offset would resolve to 0, 1, 2, 3
            sSampleIndex[s] = 
                floorInt(center - filterWidth + 0.5f);
            float weightSum = 0.0f;
            int wOffset = s * nSamples;
            for (int i = 0; i < nSamples; ++i) {
                float p = sSampleIndex[s] + 0.5f + i;
                sSampleWeight[wOffset + i] = gaussian(p - center, filterWidth);
                weightSum += sSampleWeight[wOffset + i];
            }
            float invW = 1.0f / weightSum;
            for (int i = 0; i < nSamples; ++i) {
                sSampleWeight[wOffset + i] *= invW;
            }
        }

        float* tSampleWeight = new float[dstHeight * nSamples];
        int* tSampleIndex = new int[dstHeight];
        for (int t = 0; t < dstHeight; ++t) {
            float center = ((float)t + 0.5f) / dstHeight * srcHeight;
            tSampleIndex[t] = 
                floorInt(center - filterWidth + 0.5f);
            float weightSum = 0.0f;
            int wOffset = t * nSamples;
            for (int i = 0; i < nSamples; ++i) {
                float p = tSampleIndex[t] + 0.5f + i;
                tSampleWeight[wOffset + i] = gaussian(p - center, filterWidth);
                weightSum += tSampleWeight[wOffset + i];
            }
            float invW = 1.0f / weightSum;
            for (int i = 0; i < nSamples; ++i) {
                tSampleWeight[wOffset + i] *= invW;
            }
        }

        T* dstBuffer = new T[dstWidth * dstHeight];
        for (int t = 0; t < dstHeight; ++t) {
            for (int s = 0; s < dstWidth; ++s) {
                int index = t * dstWidth + s;
                dstBuffer[index] = T(0.0f);
                for (int i = 0; i < nSamples; ++i) {
                    int srcT = clamp(tSampleIndex[t] + i, 0,
                        srcHeight - 1);
                    for (int j = 0; j < nSamples; ++j) {
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
        string type = params.getString("mapping", "uv");
        if (type == "uv") {
            Vector2 scale = params.getVector2("scale", Vector2(1.0f, 1.0f));
            Vector2 offset = params.getVector2("offset", Vector2::Zero);
            m = new UVMapping(scale, offset);
        } else if (type == "spherical") {
            Transform toTex = getTransform(params);
            m = new SphericalMapping(toTex);
        } else {
            cerr << "undefined mapping type " << type << endl;
            m = new UVMapping(Vector2(1.0f, 1.0f), Vector2::Zero);
        }
        return m;
    }

    Texture<float>* FloatConstantTextureCreator::create(const ParamSet& params,
        const SceneCache& sceneCache) const {
        float f = params.getFloat("float", 0.5f);
        return new ConstantTexture<float>(f);
    }

    Texture<float>* FloatCheckboardTextureCreator::create(
        const ParamSet& params, const SceneCache& sceneCache) const {
        TextureMapping* mapping = getTextureMapping(params);
        FloatTexturePtr T1 = sceneCache.getFloatTexture(
            params.getString("texture1"));
        FloatTexturePtr T2 = sceneCache.getFloatTexture(
            params.getString("texture2"));
        bool filter = params.getBool("filter", false);
        return new CheckboardTexture<float>(mapping, T1, T2, filter);
    }

    Texture<float>* FloatScaleTextureCreator::create(const ParamSet& params,
        const SceneCache& sceneCache) const {
        string textureName = params.getString("texture");
        string scaleName = params.getString("scale");
        FloatTexturePtr s = sceneCache.getFloatTexture(scaleName);
        FloatTexturePtr t = sceneCache.getFloatTexture(textureName);
        return new ScaleTexture<float>(t, s);
    }

    struct ImageTextureParams {
        string filePath;
        TextureMapping* mapping;
        FilterType filter;
        AddressMode addressMode;
        float gamma;
        ImageChannel channel;
        float maxAnisotropy;
    };

    static ImageTextureParams getImageTextureParams(const ParamSet& params,
        const SceneCache& sceneCache) {
        ImageTextureParams ip;
        ip.mapping = getTextureMapping(params);

        string filename = params.getString("file");
        ip.filePath = sceneCache.resolvePath(filename);

        string filterStr = params.getString("filter", "nearest");
        if (filterStr == "nearest") {
            ip.filter = FilterNone;
        } else if (filterStr == "bilinear") {
            ip.filter = FilterBilinear;
        } else if (filterStr == "trilinear") {
            ip.filter = FilterTrilinear;
        } else if (filterStr == "EWA") {
            ip.filter = FilterEWA;
        } else {
            cerr << "unrecognize filter: " << filterStr << endl;
            ip.filter = FilterNone;
        }

        string addressStr = params.getString("address", "repeat");
        if (addressStr == "repeat") {
            ip.addressMode = AddressRepeat;
        } else if (addressStr == "clamp") {
            ip.addressMode = AddressClamp;
        } else if (addressStr == "border") {
            ip.addressMode = AddressBorder;
        } else {
            cerr << "unrecognize address mode: " << addressStr << endl;
            ip.addressMode = AddressRepeat;
        }

        ip.gamma = params.getFloat("gamma", 1.0f);

        string channelStr = params.getString("channel", "All");
        if (channelStr == "R") {
            ip.channel = ChannelR;
        } else if (channelStr == "G") {
            ip.channel = ChannelG;
        } else if (channelStr == "B") {
            ip.channel = ChannelB;
        } else if (channelStr == "A") {
            ip.channel = ChannelA;
        } else if (channelStr == "All") {
            ip.channel = ChannelAll;
        } else {
            cerr << "unrecognize channel: " << channelStr << endl;
            ip.channel = ChannelAll;
        }
        ip.maxAnisotropy = params.getFloat("max_anisotropy", 10.0f);
        return ip;
    }

    Texture<float>* FloatImageTextureCreator::create(const ParamSet& params,
        const SceneCache& sceneCache) const {
        ImageTextureParams ip = getImageTextureParams(params, sceneCache);
        return new ImageTexture<float>(ip.filePath, ip.mapping, ip.filter,
            ip.addressMode, ip.gamma, ip.channel, ip.maxAnisotropy);
    }

    Texture<Color>* ColorConstantTextureCreator::create(const ParamSet& params,
        const SceneCache& sceneCache) const {
        Vector3 c = params.getVector3("color");
        return new ConstantTexture<Color>(Color(c[0], c[1], c[2]));
    }

    Texture<Color>* ColorCheckboardTextureCreator::create(
        const ParamSet& params, const SceneCache& sceneCache) const {
        TextureMapping* mapping = getTextureMapping(params);
        ColorTexturePtr T1 = sceneCache.getColorTexture(
            params.getString("texture1"));
        ColorTexturePtr T2 = sceneCache.getColorTexture(
            params.getString("texture2"));
        bool filter = params.getBool("filter", false);
        return new CheckboardTexture<Color>(mapping, T1, T2, filter);
    }

    Texture<Color>* ColorScaleTextureCreator::create(const ParamSet& params,
        const SceneCache& sceneCache) const {
        string textureName = params.getString("texture");
        string scaleName = params.getString("scale");
        FloatTexturePtr scale = sceneCache.getFloatTexture(scaleName);
        ColorTexturePtr texture = sceneCache.getColorTexture(textureName);
        return new ScaleTexture<Color>(texture, scale);
    }

    Texture<Color>* ColorImageTextureCreator::create(const ParamSet& params,
        const SceneCache& sceneCache) const {
        ImageTextureParams ip = getImageTextureParams(params, sceneCache);
        return new ImageTexture<Color>(ip.filePath, ip.mapping, ip.filter,
            ip.addressMode, ip.gamma, ip.channel);
    }

    template struct ImageBuffer<float>;
    template struct ImageBuffer<Color>; 
    template class MIPMap<float>;
    template class MIPMap<Color>; 
    template class Texture<float>;
    template class Texture<Color>;
    template class ConstantTexture<float>;
    template class ConstantTexture<Color>;
    template class CheckboardTexture<float>;
    template class CheckboardTexture<Color>;
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

}
