#include "GoblinObjMesh.h"
#include "GoblinTriangle.h"
#include "GoblinScene.h"
#include "GoblinParamSet.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

namespace Goblin {

struct TriIndex {
    int vertex;
    int normal;
    int texCoord;
    bool operator==(const TriIndex& rhs) const {
        return vertex == rhs.vertex &&
            normal == rhs.normal &&
            texCoord == rhs.texCoord;
    }

    bool operator<(const TriIndex& rhs) const {
        if (vertex == rhs.vertex) {
            if (normal == rhs.normal) {
                return texCoord < rhs.texCoord;
            } else {
                return normal < rhs.normal;
            }
        } else {
            return vertex < rhs.vertex;
        }
    }
};

struct Face {
    TriIndex index[3];
};

enum ObjFormat {
    VERTEX_ONLY = 1 << 0,
    VERTEX_UV = 1 << 1,
    VERTEX_NORMAL = 1 << 2,
    VERTEX_UV_NORMAL = 1 << 3
};

ObjMesh::ObjMesh(const std::string& filename) :
    mFilename(filename), mArea(0.0f), 
    mHasNormal(false), mHasTexCoord(false) {}

ObjMesh::~ObjMesh() {
}

void ObjMesh::init() {
    geometryCache[getId()] = this;
    load();
}

bool ObjMesh::load() {
    std::ifstream file(mFilename.c_str());
    if (!file.is_open()) {
        std::cerr << "Error can't open obj file: " 
            << mFilename << " for mesh loading" << std::endl;
        return false;
    }

    typedef std::vector<Vector3> VertexList;
    typedef std::vector<Vector3> NormalList;
    typedef std::vector<Vector2> UVList;
    typedef std::vector<Face> FaceList;
    typedef std::map<TriIndex, unsigned int> VertexMap;
    VertexMap vMap;

    VertexList vertexList;
    NormalList normalList;
    UVList uvList;
    FaceList faceList;

    ObjFormat format = VERTEX_ONLY;
    mHasNormal = false;
    mHasTexCoord = false;
    std::string line;
    int lineNum = 0;
    while(std::getline(file, line)) {
        lineNum++;
        std::stringstream stream(line);
        std::string token;
        stream >> token;
        if (token == "v") {
            Vector3 position;
            stream >> position.x >> position.y >> position.z;
            if (stream.fail()) {
                std::cout << line << std::endl;
                std::cerr << "position syntax error on line " 
                    << lineNum << std::endl;
                return false;
            }
            vertexList.push_back(position);
        } else if (token == "vn") {
            Vector3 normal;
            stream >> normal.x >> normal.y >> normal.z;
            if (stream.fail()) {
                std::cout << line << std::endl;
                std::cerr << "normal syntax error on line " 
                    << lineNum << std::endl;
                return false;
            }
            normalList.push_back(normal);
        } else if (token == "vt") {
            Vector2 uv;
            stream >> uv.x >> uv.y;
            if (stream.fail()) {
                std::cerr << "uv syntax error on line " 
                    <<lineNum << std::endl;
                return false;
            }
            uvList.push_back(uv);
        } else if (token == "f") {
            std::string faceToken;
            std::vector<std::string> faceTokens;
            while(true) {
                stream >> faceToken;
                if (stream.fail()) {
                    break;
                }
                faceTokens.push_back(faceToken);
            }

            size_t faceTokensNum = faceTokens.size();
            if (faceTokensNum > 4 || faceTokensNum < 3) {
                std::cerr << "incorrect face vertices number on line " 
                    << lineNum << std::endl;
                return false;
            }

            // the first face read in, figure out the scan format
            if (faceList.size() == 0) {
                std::string faceToken = faceTokens[0];
                if (faceToken.find("//") != std::string::npos) {
                    format = VERTEX_NORMAL;
                    mHasNormal = true;
                } else if (faceToken.find('/') == std::string::npos) {
                    format = VERTEX_ONLY;
                } else {
                    size_t p1 = faceToken.find('/');
                    size_t p2 = faceToken.rfind('/');
                    if (p1 == p2) {
                        format = VERTEX_UV;
                        mHasTexCoord = true;
                    } else {
                        format = VERTEX_UV_NORMAL;
                        mHasNormal = true;
                        mHasTexCoord = true;
                    }
                }
            }
            // now scan in the face
            TriIndex triIndex[4];
            for (size_t i = 0; i < faceTokensNum; ++i) {
                const char* token = faceTokens[i].c_str();
                switch(format) {
                case VERTEX_ONLY:
                    triIndex[i].vertex = atoi(token);
                    triIndex[i].normal = 0;
                    triIndex[i].texCoord = 0;
                    break;
                case VERTEX_UV:
                    triIndex[i].vertex = atoi(token);
                    token += strcspn(token, "/") + 1;
                    triIndex[i].texCoord = atoi(token);
                    triIndex[i].normal = 0;
                    break;
                case VERTEX_NORMAL:
                    triIndex[i].vertex = atoi(token);
                    token += strcspn(token, "/") + 2;
                    triIndex[i].normal = atoi(token);
                    triIndex[i].texCoord = 0;
                    break;
                case VERTEX_UV_NORMAL:
                    triIndex[i].vertex = atoi(token);
                    token += strcspn(token, "/") + 1;
                    triIndex[i].texCoord = atoi(token);
                    token += strcspn(token, "/") + 1;
                    triIndex[i].normal = atoi(token);
                    break;
                default:
                    std::cerr << "unrecognize face format on line " 
                        << lineNum << std::endl;
                    break;
                }
            }

            Face f1 = {{triIndex[0], triIndex[1], triIndex[2]}};
            faceList.push_back(f1);
            if (faceTokensNum == 4) {
                Face f2 = {{triIndex[0], triIndex[2], triIndex[3]}};
                faceList.push_back(f2);
            }
        }
        //std::cout <<line << std::endl;
        line.clear();
    }

    int verticesNum = static_cast<int>(vertexList.size());
    int normalNum = static_cast<int>(normalList.size());
    int texCoordNum = static_cast<int>(uvList.size());
        
    for (size_t i = 0; i < faceList.size(); ++i) {
        Face& face = faceList[i];
        for (size_t j = 0; j < 3; ++j) {
            int& vIndex = face.index[j].vertex;
            int& nIndex = face.index[j].normal;
            int& tIndex = face.index[j].texCoord;
            // some obj do backward index like what python do
            if (vIndex < 0) {
                vIndex += verticesNum + 1;
            }
            if (nIndex < 0) {
                nIndex += normalNum + 1;
            }
            if (tIndex < 0) {
                tIndex += texCoordNum + 1;
            }
            // obj face indexing starts on 1 instead of 0
            vIndex--;
            nIndex--;
            tIndex--;
            if (vIndex < 0 || vIndex >= verticesNum ||
                nIndex < -1 || nIndex >= normalNum ||
                tIndex < -1 || tIndex >= texCoordNum) {
                    std::cerr << "invalid index in face " << i 
                        << std::endl;
                    return false;
            }
        }
    }
    //mVertices.reserve(faceList.size());
    mTriangles.reserve(faceList.size());

    unsigned int vIndexCounter = 0;
    for (size_t i = 0; i < faceList.size(); ++i) {
        const Face& face = faceList[i];
        TriangleIndex triangle;
        for (size_t j = 0; j < 3; ++j) {
            std::pair<VertexMap::iterator, bool> rv = 
                vMap.insert(std::make_pair(face.index[j], vIndexCounter)); 
            // new inserted VertexIndex
            if (rv.second) {
                Vertex v;
                int vIndex = face.index[j].vertex;
                v.position = vertexList[vIndex];
                int nIndex = face.index[j].normal;
                v.normal = nIndex == -1 ? Vector3::Zero : normalList[nIndex];
                int tIndex = face.index[j].texCoord;
                v.texC = tIndex == -1 ? Vector2::Zero : uvList[tIndex];
                mVertices.push_back(v);
                mBBox.expand(v.position);
                vIndexCounter++;
            }
            // (vIndex, nIndex, tIndex) map to the final triangle index
            triangle.v[j] = rv.first->second;
        }
        mTriangles.push_back(triangle);
    }
    recalculateArea();
    std::cout << "Successfully loaded mesh '" << mFilename << "'.\n";
    return true;
}

bool ObjMesh::intersect(const Ray& ray) const {
    return false;
}

bool ObjMesh::intersect(const Ray& ray, float* epsilon, 
    Fragment* fragment) const {
    return false;
}

BBox ObjMesh::getObjectBound() const {
    return mBBox;
}

void ObjMesh::refine(GeometryList& refinedGeometries) const {
    size_t faceNum = mTriangles.size();
    if (mRefinedMeshes.size() != faceNum) {
        mRefinedMeshes.clear();
        mRefinedMeshes.resize(faceNum, Triangle(this));
        for (size_t i = 0; i < faceNum; ++i) {
            mRefinedMeshes[i].setIndex(i);
        }
    }
    for (size_t i = 0; i < faceNum; ++i) {
        refinedGeometries.push_back(&mRefinedMeshes[i]);
    }
}

void ObjMesh::recalculateArea() {
    mArea = 0.0f;
    for (size_t i = 0; i < mTriangles.size(); ++i) {
        size_t i0 = mTriangles[i].v[0];
        size_t i1 = mTriangles[i].v[1];
        size_t i2 = mTriangles[i].v[2];
        Vector3 v0 = mVertices[i0].position;
        Vector3 v1 = mVertices[i1].position;
        Vector3 v2 = mVertices[i2].position;
        mArea += 0.5f * length(cross(v1 - v0, v2 - v0));
    }
}

Geometry* createPolygonMesh(const ParamSet& params, const SceneCache& sceneCache) {
	std::string filename = params.getString("file");
	std::string filePath = sceneCache.resolvePath(filename);
    return new ObjMesh(filePath);
}

}
