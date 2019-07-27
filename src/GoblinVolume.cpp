#include "GoblinVolume.h"
#include "GoblinParamSet.h"
#include "GoblinRay.h"
#include "GoblinScene.h"

#include <fstream>
#include <sstream>
#include <cstring>

namespace Goblin {

bool VolumeRegion::intersect(const Ray& ray, 
    float* tMin, float* tMax) const {
    return mLocalRegion.intersect(mToWorld.invertRay(ray), tMin, tMax);
}

float VolumeRegion::phase(const Vector3& p, const Vector3& wi, 
    const Vector3& wo) const {
    if (!mLocalRegion.contain(mToWorld.invertPoint(p))) {
        return 0.0f;
    }
    return phaseHG(wi, wo, mG);
}

Color HomogeneousVolumeRegion::transmittance(const Ray& ray,
    const RNG& rng) const {
    // homogeneous transmittance can be analytically solved by Beer's law
    float tMin, tMax;
    if (!mLocalRegion.intersect(mToWorld.invertRay(ray), &tMin, &tMax)) {
        return Color(1.0f);
    }
    Color tau = length(ray(tMax) - ray(tMin)) * mAttenuation;;
    return Color(exp(-tau.r), exp(-tau.g), exp(-tau.b));
}

enum VolumeEncoding {
    VolumeEncodingFloat32 = 1,
    VolumeEncodingFloat16 = 2,
    VolumeEncodingUInt8 = 3,
    VolumeEncodingQuantizedDirections = 4
};

// The .vol file is Mitsuba's binary exchange format storing 3d grid data
// See Mitsuba documentation for the detail binary layout spec
// Bytes 1-3    ASCII Bytes 'V', 'O', and 'L'
// Byte  4      File format version number (currently 3)
// Bytes 5-8    Encoding identifier (32-bit integer).
//              The following choices are available:
//              1. Dense float32-based representation
//              2. Dense float16-based representation
//              3. Dense uint8-based representation
//                 (The range 0..255 will be mapped to 0..1)
//              4. Dense quantized directions. The directions are stored
//                 in spherical coordinates with a total storage cost of
//                 16 bit per entry.
// Bytes 9-12   Number of cells along the X axis (32 bit integer)
// Bytes 13-16  Number of cells along the Y axis (32 bit integer)
// Bytes 17-20  Number of cells along the Z axis (32 bit integer)
// Bytes 21-24  Number of channels (32 bit integer, supported values: 1 or 3)
// Bytes 25-48  Axis-aligned bounding box of the data stored in
//              single precision (order: xmin, ymin, zmin, xmax, ymax, zmax)
// Bytes 49-*   Binary data of the volume stored in the specified encoding.
//              The data are ordered so that the following C-style indexing
//              operation makes sense after the file has been
//              mapped into memory:
//              data[((zpos*yres + ypos)*xres + xpos)*channels + chan]
//              where (xpos, ypos, zpos, chan) denotes the lookup location.
struct VolHeader
{
    bool isValid(std::string* msg = nullptr) const {
        if (mSignature[0] != 'V' ||
            mSignature[1] != 'O' ||
            mSignature[2] != 'L') {
            if (msg) {
                *msg = "invalid file format";
            }
            return false;
        }
        if (mNChannel != 1 && mNChannel != 3) {
            if (msg) {
                *msg = "invalid channel count (only support 1 or 3)";
            }
            return false;
        }
        if (mNx <= 0 || mNy <= 0 || mNz <= 0) {
            if (msg) {
                *msg = "invalid grid resolution";
            }
            return false;
        }
        if (mEncoding < VolumeEncodingFloat32 ||
            mEncoding > VolumeEncodingQuantizedDirections) {
            if (msg){
                *msg = "invalid encoding";
            }
            return false;
        }
        return true;
    }

    void print() const {
        std::string encodingStr;
        if (mEncoding == VolumeEncodingFloat32) {
            encodingStr = "Float32";
        } else if (mEncoding == VolumeEncodingFloat16) {
            encodingStr = "Float16";
        } else if (mEncoding == VolumeEncodingUInt8) {
            encodingStr = "UInt8";
        } else if (mEncoding == VolumeEncodingQuantizedDirections) {
            encodingStr = "Quantized Directions";
        } else {
            encodingStr = "Invalid encoding";
        }
        std::cout << "encoding: " <<encodingStr << "\n" <<
            "nX: " << mNx << " nY " << mNy << " nZ " << mNz << "\n" <<
            "bbox min: " << mMin << "\n" <<
            "bbox max: " << mMax << std::endl;
    }


    size_t getEntryCount() const {
        return mNx * mNy * mNz;
    }

    char mSignature[4];
    int32_t mEncoding;
    int32_t mNx;
    int32_t mNy;
    int32_t mNz;
    int32_t mNChannel;
    Vector3 mMin;
    Vector3 mMax;
};

class VolumeGrid
{
public:
    VolumeGrid(int32_t nx, int32_t ny, int32_t nz,
        int32_t nChannel, const BBox& bbox, float* data):
        mNx(nx), mNy(ny), mNz(nz), mNChannel(nChannel), mBBox(bbox),
        mVoxelData(nx * ny * nz * nChannel) {
        Vector3 dim = mBBox.pMax - mBBox.pMin;
        mNormalizeTerm = Vector3(1.0f / dim.x, 1.0f / dim.y, 1.0f / dim.z);
        memcpy(&mVoxelData[0], data, mVoxelData.size() * sizeof(float));
    }

    Color getVoxel(int x, int y, int z) const {
        if (x < 0 || x >= mNx ||
            y < 0 || y >= mNy ||
            z < 0 || z >= mNz) {
            return Color(0.0f);
        } else if (mNChannel == 1) {
            int offset = z * mNx * mNy +y * mNx + x;
            return Color(mVoxelData[offset]);
        } else if (mNChannel == 3) {
            int offset = 3 * (z * mNx * mNy +y * mNx + x);
            return Color(mVoxelData[offset],
                mVoxelData[offset + 1],
                mVoxelData[offset + 2]);
        } else {
            return Color(0.0f);
        }
    }

    Color eval(const Vector3& pLocal) const {
        Vector3 fIndex = (pLocal - mBBox.pMin);
        // mapping the position inside bounding box to
        // (0, 0, 0) - (nx -1, ny - 1, nz - 1)
        fIndex.x *= mNormalizeTerm.x * mNx - 0.5f;
        fIndex.y *= mNormalizeTerm.y * mNy - 0.5f;
        fIndex.z *= mNormalizeTerm.z * mNz - 0.5f;
        // trilinear interpolate the voxel lookup results
        int ix = floorInt(fIndex.x);
        int iy = floorInt(fIndex.y);
        int iz = floorInt(fIndex.z);
        float dx = fIndex.x - ix;
        float dy = fIndex.y - iy;
        float dz = fIndex.z - iz;
        Color d00 = lerp(dx,
            getVoxel(ix, iy, iz), getVoxel(ix + 1, iy, iz));
        Color d10 = lerp(dx,
            getVoxel(ix, iy + 1, iz), getVoxel(ix + 1, iy + 1, iz));
        Color d01 = lerp(dx,
            getVoxel(ix, iy, iz + 1), getVoxel(ix + 1, iy, iz + 1));
        Color d11 = lerp(dx,
            getVoxel(ix, iy + 1, iz + 1), getVoxel(ix + 1, iy + 1, iz + 1));
        Color d0 = lerp(dy, d00, d10);
        Color d1 = lerp(dy, d01, d11);
        return lerp(dz, d0, d1);
    }

    const BBox& getBBox() const { return mBBox; }

    void fillHeaderInfo(VolHeader& header) const {
        header.mSignature[0] = 'V';
        header.mSignature[1] = 'O';
        header.mSignature[2] = 'L';
        header.mSignature[3] = 3;
        header.mEncoding = VolumeEncodingFloat32;
        header.mNx = mNx;
        header.mNy = mNy;
        header.mNz = mNz;
        header.mNChannel = mNChannel;
        header.mMin = mBBox.pMin;
        header.mMax = mBBox.pMax;
    }

    const float* getData() const {
        return &mVoxelData[0];
    }

private:
    int32_t mNx;
    int32_t mNy;
    int32_t mNz;
    int32_t mNChannel;
    BBox mBBox;
    Vector3 mNormalizeTerm;
    std::vector<float> mVoxelData;
};

VolumeGrid* loadVolFile(const std::string& filePath, std::string* error = nullptr) {
    VolumeGrid* result = nullptr;
    std::ifstream stream(filePath.c_str(), std::ios::in | std::ios::binary);
    if (stream.is_open()) {
        std::streampos begin = stream.tellg();
        stream.seekg (0, std::ios::end);
        std::streampos end = stream.tellg();
        std::cout << "successfully open vol file " << filePath << std::endl;
        stream.seekg(0);
        VolHeader header;
        stream.read(reinterpret_cast<char*>(&header), sizeof(header));
        std::string msg;
        if (header.isValid(&msg)) {
            header.print();
            // TODO suppport loading other encoding type
            if (header.mEncoding == VolumeEncodingFloat32) {
                size_t nFloat = header.getEntryCount() * header.mNChannel;
                stream.seekg(sizeof(header));
                float* data = new float[nFloat];
                stream.read(reinterpret_cast<char*>(data),
                    nFloat * sizeof(float));
                result = new VolumeGrid(header.mNx, header.mNy, header.mNz,
                    header.mNChannel, BBox(header.mMin, header.mMax), data);
                delete[] data;
            }
        } else {
            std::cerr << msg << std::endl;
        }
    } else {
        std::cout << "fail to open vol file " << filePath <<std::endl;
    }
    stream.close();
    return result;
}

bool writeVolFile(const std::string& filePath,
    const VolumeGrid* volumeGrid, std::string* error = nullptr) {
    bool result = false;
    std::ofstream stream(filePath.c_str(), std::ios::out | std::ios::binary);
    if (stream.is_open()) {
        stream.seekp(0);
        VolHeader header;
        volumeGrid->fillHeaderInfo(header);
        stream.write(reinterpret_cast<char*>(&header), sizeof(header));
        std::string msg;
        size_t nFloat = header.getEntryCount() * header.mNChannel;
        stream.seekp(sizeof(header));
        stream.write(reinterpret_cast<const char*>(volumeGrid->getData()),
            nFloat * sizeof(float));
        result = true;
    } else {
        if (error != nullptr){
            *error = "fail to open vol file ";
        }
    }
    stream.close();
    return result;
}

HeterogeneousVolumeRegion::HeterogeneousVolumeRegion(
    VolumeGrid* density, const Color& albedo, float g,
    float stepSize, int sampleNum, const BBox& b, const Transform& toWorld):
    VolumeRegion(g, stepSize, sampleNum, b, toWorld, false),
    mDensity(density), mAlbedo(albedo)
{}

HeterogeneousVolumeRegion::~HeterogeneousVolumeRegion() {
    if (mDensity) {
        delete mDensity;
        mDensity = nullptr;
    }
}

void HeterogeneousVolumeRegion::eval(const Vector3& p,
    Color& attenuation, Color& scatter, Color& emission) const {
    Vector3 pLocal = mToWorld.invertPoint(p);
    bool inside = mLocalRegion.contain(pLocal);
    if (inside) {
        attenuation = (mDensity != nullptr) ?
            mDensity->eval(pLocal) : Color(0.0f);
        scatter = attenuation * mAlbedo;
        emission = Color(0.0f);
    } else {
        attenuation = Color(0.0f);
        scatter = Color(0.0f);
        emission = Color(0.0f);
    }
}

Color HeterogeneousVolumeRegion::getAttenuation(const Vector3& p) const {
    Vector3 pLocal = mToWorld.invertPoint(p);
    bool inside = mLocalRegion.contain(pLocal);
    if (inside) {
        return (mDensity != nullptr) ? mDensity->eval(pLocal) : Color(0.0f);
    } else {
        return Color(0.0f);
    }
}

Color HeterogeneousVolumeRegion::transmittance(const Ray& ray,
    const RNG& rng) const {
    // TODO replace biased ray marching approach with unbiased solution
    // (ratio tracking or residual ratio tracking)
    float tMin, tMax;
    if (!mLocalRegion.intersect(mToWorld.invertRay(ray), &tMin, &tMax)) {
        return Color(0.0f);
    }
    float stepSize = getSampleStepSize();
    float tCurrent = tMin;
    float jitter = rng.randomFloat() * stepSize;
    Color tau = jitter * getAttenuation(ray(tCurrent));
    tCurrent += jitter;
    while (tCurrent + stepSize < tMax) {
        tau += stepSize * getAttenuation(ray(tCurrent));
        tCurrent += stepSize;
    }
    tau += (tMax - tCurrent) * getAttenuation(ray(tCurrent));
    return Color(exp(-tau.r), exp(-tau.g), exp(-tau.b));
}

VolumeRegion* createHomogeneousVolume(
    const ParamSet& params, const SceneCache& sceneCache) {
    Vector3 attenuation = params.getVector3("attenuation");
    Vector3 albedo = params.getVector3("albedo");
    Vector3 emission = params.getVector3("emission");
    float g = params.getFloat("g", 0.0f);
    int sampleNum = params.getInt("sample_num", 5);
    Vector3 vMin = params.getVector3("box_min");
    Vector3 vMax = params.getVector3("box_max");
    BBox b(vMin, vMax);
    Transform toWorld = getTransform(params);
    return new HomogeneousVolumeRegion(
		Color(attenuation[0], attenuation[1], attenuation[2]),
		Color(albedo[0], albedo[1], albedo[2]),
		Color(emission[0], emission[1], emission[2]),
        g, sampleNum, b, toWorld);
}

VolumeRegion* createHeterogeneousVolume(
    const ParamSet& params, const SceneCache& sceneCache) {
    std::string densityPath = sceneCache.resolvePath(
        params.getString("density_grid"));
    std::string error;
    VolumeGrid* densityGrid = loadVolFile(densityPath, &error);
    if (densityGrid == nullptr) {
        std::cerr << error << std::endl;
        // construct a one entry default VolumeGid as fallback
        float d = 1.0f;
        densityGrid = new VolumeGrid(1, 1, 1, 1,
            BBox(Vector3(-1, -1, -1), Vector3(1, 1, 1)), &d);
    }

    Vector3 albedo = params.getVector3("albedo");
    float g = params.getFloat("g", 0.0f);
    float stepSize = params.getFloat("step_size", 0.1f);
    int sampleNum = params.getInt("sample_num", 5);
    BBox b = densityGrid->getBBox();
    Transform toWorld = getTransform(params);
    return new HeterogeneousVolumeRegion(densityGrid,
		Color(albedo[0], albedo[1], albedo[2]), g,
        stepSize, sampleNum, b, toWorld);
}
}
