#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <map>

#include "GoblinObjMesh.h"
namespace Goblin {

    struct TriIndex {
        int vertex;
        int normal;
        int texCoord;

        bool operator<(const TriIndex& rhs) const {
            if(vertex == rhs.vertex) {
                if(normal == rhs.normal) {
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

    bool ObjMesh::load2(const std::string& filename) {
        std::ifstream file = std::ifstream(filename.c_str());
        if(!file.is_open()) {
            std::cerr << "Error can't open obj file: " 
                << filename << " for mesh loading" << std::endl;
            return false;
        }

        std::string line;
        std::string token;
        int lineNum = 0;
        while(std::getline(file, line)) {
            lineNum++;
        }
        std::cout << "done with reading, start assembling" << std::endl;
        return true;
    }

    bool ObjMesh::load(const std::string& filename) {
        std::ifstream file = std::ifstream(filename.c_str());
        if(!file.is_open()) {
            std::cerr << "Error can't open obj file: " 
                << filename << " for mesh loading" << std::endl;
            return false;
        }

        typedef std::vector<Vector3> VertexList;
        typedef std::vector<Vector3> NormalList;
        typedef std::vector<Vector2> UVList;
        typedef std::vector<Face> FaceList;

        static const char* scanVertex = "%d";
        static const char* scanVertexUV = "%d/%d";
        static const char* scanVertexNormal = "%d//%d";
        static const char* scanVertexUVNormal = "%d/%d/%d";

        VertexList vertexList;
        NormalList normalList;
        UVList uvList;
        FaceList faceList;

        ObjFormat format = VERTEX_ONLY;
        bool hasNormal = false;
        bool hasTexCoord = false;

        std::string line;
        std::string token;
        int lineNum = 0;
        while(std::getline(file, line)) {
            lineNum++;
            std::stringstream stream(line);
            stream >> token;
            if(token == "v") {
                Vector3 position;
                stream >> position.x >> position.y >> position.z;
                if(stream.fail()) {
                    std::cerr << "position syntax error on line " 
                        << lineNum << std::endl;
                    return false;
                }
                vertexList.push_back(position);
            } else if(token == "vn") {
                Vector3 normal;
                stream >> normal.x >> normal.y >> normal.z;
                if(stream.fail()) {
                    std::cerr << "normal syntax error on line " 
                        << lineNum << std::endl;
                    return false;
                }
                normalList.push_back(normal);
            } else if(token == "vt") {
                Vector2 uv;
                stream >> uv.x >> uv.y;
                if(stream.fail()) {
                    std::cerr << "uv syntax error on line " 
                        <<lineNum << std::endl;
                    return false;
                }
                uvList.push_back(uv);
            } else if(token == "f") {
                std::string faceToken;
                std::vector<std::string> faceTokens;
                while(true) {
                    stream >> faceToken;
                    if(stream.fail()) {
                        break;
                    }
                    faceTokens.push_back(faceToken);
                }

                size_t faceTokensNum = faceTokens.size();
                if(faceTokensNum > 4 || faceTokensNum < 3) {
                    std::cerr << "incorrect face vertices number on line " 
                        << lineNum << std::endl;
                    return false;
                }

                // the first face read in, figure out the scan format
                if(faceList.size() == 0) {
                    std::string faceToken = faceTokens[0];
                    if(faceToken.find("//") != std::string::npos) {
                        format = VERTEX_NORMAL;
                        hasNormal = true;
                    } else if(faceToken.find('/') == std::string::npos) {
                        format = VERTEX_ONLY;
                    } else {
                        size_t p1 = faceToken.find('/');
                        size_t p2 = faceToken.rfind('/');
                        if(p1 == p2) {
                            format = VERTEX_UV;
                            hasTexCoord = true;
                        } else {
                            format = VERTEX_UV_NORMAL;
                            hasNormal = true;
                            hasTexCoord = true;
                        }
                    }
                }
                // now scan in the face
                TriIndex triIndex[4];
                for(size_t i = 0; i < faceTokensNum; ++i) {
                    switch(format) {
                    case VERTEX_ONLY:
                        sscanf(faceTokens[i].c_str(), scanVertex, 
                            &triIndex[i].vertex);
                        triIndex[i].normal = 0;
                        triIndex[i].texCoord = 0;
                        break;
                    case VERTEX_UV:
                        sscanf(faceTokens[i].c_str(), scanVertexUV,
                            &triIndex[i].vertex, &triIndex[i].texCoord);
                        triIndex[i].normal = 0;
                        break;
                    case VERTEX_NORMAL:
                        sscanf(faceTokens[i].c_str(), scanVertexNormal,
                            &triIndex[i].vertex, &triIndex[i].normal);
                        triIndex[i].texCoord = 0;
                        break;
                    case VERTEX_UV_NORMAL:
                        sscanf(faceTokens[i].c_str(), scanVertexUVNormal,
                            &triIndex[i].vertex, &triIndex[i].texCoord,
                            &triIndex[i].normal);
                        break;
                    default:
                        std::cerr << "unrecognize face format on line " 
                            << lineNum << std::endl;
                        break;
                    }
                }
                // obj face indexing starts on 1 instead of 0
                for(size_t i = 0; i < faceTokensNum; ++i) {
                    triIndex[i].vertex--;
                    triIndex[i].normal--;
                    triIndex[i].texCoord--;
                }
                Face f1 = {{triIndex[0], triIndex[1], triIndex[2]}};
                faceList.push_back(f1);
                if(faceTokensNum == 4) {
                    Face f2 = {{triIndex[0], triIndex[2], triIndex[3]}};
                    faceList.push_back(f2);
                }
            }
            //std::cout <<line << std::endl;
            line.clear();
        }

        int verticesNum = vertexList.size();
        int normalNum = normalList.size();
        int texCoordNum = uvList.size();
        
        for(size_t i = 0; i < faceList.size(); ++i) {
            const Face& face = faceList[i];
            for(size_t j = 0; j < 3; ++j) {
                int vIndex = face.index[j].vertex;
                int nIndex = face.index[j].normal;
                int tIndex = face.index[j].texCoord;
                if(vIndex < 0 || vIndex >= verticesNum ||
                    nIndex < -1 || nIndex >= normalNum ||
                    tIndex < -1 || tIndex >= texCoordNum) {
                        std::cerr << "invalid index in face " << i 
                            << std::endl;
                        return false;
                }
            }
        }
        std::cout << "done with reading, start assembling" << std::endl;
        vertices.clear();
        vertices.reserve(vertexList.size());
        triangles.clear();
        triangles.reserve(faceList.size());
        typedef std::map<TriIndex, unsigned int> VertexMap;
        VertexMap vMap;
        unsigned int vIndexCounter = 0;
        for(size_t i = 0; i < faceList.size(); ++i) {
            const Face& face = faceList[i];
            MeshTriangle triangle;
            for(size_t j = 0; j < 3; ++j) {
                std::pair<VertexMap::iterator, bool> rv = 
                    vMap.insert(std::make_pair(face.index[j], vIndexCounter)); 
                // new inserted VertexIndex
                if(rv.second) {
                    Vertex v;
                    int vIndex = face.index[j].vertex;
                    v.position = vertexList[vIndex];
                    int nIndex = face.index[j].normal;
                    v.normal = nIndex == -1 ? Vector3::Zero : normalList[nIndex];
                    int tIndex = face.index[j].texCoord;
                    v.texC = tIndex == -1 ? Vector2::Zero : uvList[tIndex];
                    vertices.push_back(v);
                    vIndexCounter++;
                }
                // (vIndex, nIndex, tIndex) map to the final triangle index
                triangle.v[j] = rv.first->second;
            }
            triangles.push_back(triangle);
        }

        std::cout << "Successfully loaded mesh '" << filename << "'.\n";
        return true;
    }
}