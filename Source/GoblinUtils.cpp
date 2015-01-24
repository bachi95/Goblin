#include "GoblinUtils.h"
#include "GoblinVector.h"
#include "GoblinColor.h"

#include <ctime>
#include <limits>
#include <boost/generator_iterator.hpp>
#include <boost/random.hpp>

namespace Goblin {

    typedef boost::mt19937 RNGType;
    typedef boost::uniform_real<float> RealDist;
    typedef boost::uniform_int<uint32_t> UInt32Dist;
    typedef boost::variate_generator< RNGType, boost::uniform_real<float> >
        RealGenerator;
    typedef boost::variate_generator< RNGType, boost::uniform_int<uint32_t> >
        UInt32Generator;

    class RNGImp {
    public:
        RNGImp();
        ~RNGImp();
        float randomFloat() const;
        uint32_t randomUInt() const;
    private:
        RNGType* mEngine;
        RealDist* mRealDist;
        UInt32Dist* mUInt32Dist;
        RealGenerator* mRealGenerator;
        UInt32Generator* mUInt32Generator;
    };

    RNGImp::RNGImp() {
        mEngine = new RNGType(static_cast<uint32_t>(rand()));
        mRealDist = new boost::uniform_real<float>(0.0f, 1.0f);
        mUInt32Dist = new boost::uniform_int<uint32_t>(0, 
            numeric_limits<uint32_t>::max() );
        mRealGenerator = new RealGenerator(*mEngine, *mRealDist);
        mUInt32Generator = new UInt32Generator(*mEngine, *mUInt32Dist);
    }

    RNGImp::~RNGImp() {
        delete mEngine;
        delete mRealDist;
        delete mUInt32Dist;
        delete mRealGenerator;
        delete mUInt32Generator;
    }

    float RNGImp::randomFloat() const {
        return (*mRealGenerator)();
    }

    uint32_t RNGImp::randomUInt() const {
        return (*mUInt32Generator)();
    }

    RNG::RNG() {
        mRNGImp = new RNGImp();
    }

    RNG::~RNG() {
        if(mRNGImp) {
            delete mRNGImp;
            mRNGImp = NULL;
        }
    }

    float RNG::randomFloat() const {
        return mRNGImp->randomFloat();
    }

    uint32_t RNG::randomUInt() const {
        return mRNGImp->randomUInt();
    }

    void coordinateAxises(const Vector3& a1, Vector3* a2, Vector3* a3) {
        // in case you throw in case like a1 = Vector3(0, 1, 0)
        if(fabsf(a1.x) > fabsf(a1.y)) {
            float invLen = 1.0f / sqrtf(a1.x * a1.x + a1.z * a1.z);
            // dot(a1, a2) = 0 <-> penpenticular to each other
            *a2 = Vector3(-a1.z * invLen, 0.0f, a1.x * invLen);
        } else {
            float invLen = 1.0f / sqrtf(a1.y * a1.y + a1.z * a1.z);
            *a2 = Vector3(0.0f, -a1.z * invLen, a1.y * invLen);
        }
        *a3 = cross(a1, *a2);
    }

    bool quadratic(float A, float B, float C, float* t1, float* t2) {
        float discriminant = B * B - 4.0f * A * C;
        if(discriminant < 0.0f) {
            return false;
        }
        float rootDiscrim = sqrt(discriminant);
        float q;
        // a small trick to avoid numeric error introduced from naive
        // implementation when B or -B close to rootDiscrim
        if(B < 0) {
            q = -0.5f * (B - rootDiscrim);
        } else {
            q = -0.5f * (B + rootDiscrim);
        }
        *t1 = q / A;
        *t2 = C / q;
        if(*t1 > *t2) {
            std::swap(*t1, *t2);
        }
        return true;
    }

    void drawLine(const Vector2& p0, const Vector2& p1, Color* buffer, 
        int xRes, int yRes, const Color& color) {
        int x0 = roundInt(p0.x);
        int y0 = roundInt(p0.y); 
        int x1 = roundInt(p1.x);
        int y1 = roundInt(p1.y); 

        int dx = x1 - x0;
        int dy = y1 - y0;
        int xStep = (dx > 0 ? 1 : (dx < 0 ? -1 : 0));
        int yStep = (dy > 0 ? 1 : (dy < 0 ? -1 : 0));
        if(dx < 0) {
            dx = -dx;
        }
        if(dy < 0) {
            dy = -dy;
        }
        
        int ax = 2 * dx;
        int ay = 2 * dy;
        // draw line with Bresenham's algorithm
        if(dx > dy) {
            int accumD = ay - dx;
            int y = y0;
            for(int x = x0; x != x1; x += xStep) {
                if(0 <= x && x < xRes && 0 <= y && y < yRes) {
                    int index = y * xRes + x;
                    float a = color.a;
                    buffer[index] = a * color + (1.0f - a) * buffer[index];
                    buffer[index].a = 1.0f;
                }
                if(accumD >= 0) {
                    accumD += ay - ax;
                    y += yStep;
                } else {
                    accumD += ay;
                }
            }
        } else {
            int accumD = ax - dy;
            int x = x0;
            for(int y = y0; y != y1; y += yStep) {
                if(0 <= x && x < xRes && 0 <= y && y < yRes) {
                    int index = y * xRes + x;
                    float a = color.a;
                    buffer[index] = a * color + (1.0f - a) * buffer[index];
                    buffer[index].a = 1.0f;
                }
                if(accumD >= 0) {
                    accumD += ax - ay;
                    x += xStep;
                } else {
                    accumD += ax;
                }
            }
        }
    }

    void drawPoint(const Vector2& p, Color* buffer, int xRes, int yRes,
        const Color& color, int radius) {
        if(radius <= 0) {
            return;
        }
        int x0 = roundInt(p.x);
        int y0 = roundInt(p.y);
        int squareRadius = radius * radius;
        for(int i = -radius; i <= radius; ++i) {
            int i2 = i * i;
            for(int j = -radius; j <= radius; ++j) {
                if(i2 + j * j > squareRadius) {
                    continue;
                }
                int x = x0 + j;
                int y = y0 + i;
                if(0 <= x && x < xRes && 0 <= y && y < yRes) {
                    int index = y * xRes + x;
                    float a = color.a;
                    buffer[index] = a * color + (1.0f - a) * buffer[index];
                    buffer[index].a = 1.0f;
                }
            }
        }
    }

}

