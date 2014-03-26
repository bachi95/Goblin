#include "GoblinImageIO.h"
#include "GoblinColor.h"
#include "GoblinUtils.h"

#include <iostream>
#include <cstdio>
#include <cassert>
#include <png.h>


#ifdef GOBLIN_ENABLE_EXR
#include <ImfInputFile.h>
#include <ImfRgbaFile.h>
#include <ImfOutputFile.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>
#include <half.h>
using namespace Imf;
using namespace Imath;
#endif //GOBLIN_ENABLE_EXR

namespace Goblin {

#ifdef GOBLIN_ENABLE_EXR
    static Color* loadImageEXR(const string& filename, 
        int *width, int *height) {
        half* rgba = NULL;
        Color* colorBuffer = NULL;
        try {
            InputFile file(filename.c_str());
            Box2i dw = file.header().dataWindow();
            *width = dw.max.x - dw.min.x + 1;
            *height = dw.max.y - dw.min.y + 1;

            rgba = new half[4 * (*width) * (*height)];
            FrameBuffer frameBuffer;
            frameBuffer.insert("R", Slice(
                HALF,                               // type
                (char*)(rgba),                      // base
                4 * sizeof(half),                   // xstride 
                4 * sizeof(half) * (*width),        // ystride
                1, 1,                               // x/y sampling
                0.0));                              // fill value
             frameBuffer.insert("G", Slice(
                HALF,                               // type
                (char*)(rgba) + 1 * sizeof(half),   // base
                4 * sizeof(half),                   // xstride 
                4 * sizeof(half) * (*width),        // ystride
                1, 1,                               // x/y sampling
                0.0));                              // fill value
             frameBuffer.insert("B", Slice(
                HALF,                               // type
                (char*)(rgba) + 2 * sizeof(half),   // base
                4 * sizeof(half),                   // xstride 
                4 * sizeof(half) * (*width),        // ystride
                1, 1,                               // x/y sampling
                0.0));                              // fill value
             frameBuffer.insert("A", Slice(
                HALF,                               // type
                (char*)(rgba) + 3 * sizeof(half),   // base
                4 * sizeof(half),                   // xstride 
                4 * sizeof(half) * (*width),        // ystride
                1, 1,                               // x/y sampling
                1.0));                              // fill value

            file.setFrameBuffer(frameBuffer); 
            file.readPixels(dw.min.y, dw.max.y);

            colorBuffer = new Color[(*width) * (*height)];
            for(int i = 0; i < (*width) * (*height); ++i) {
                colorBuffer[i].r = rgba[4 * i];
                colorBuffer[i].g = rgba[4 * i + 1];
                colorBuffer[i].b = rgba[4 * i + 2];
                colorBuffer[i].a = rgba[4 * i + 3];
            }
            delete [] rgba;
            return colorBuffer;

        } catch (const std::exception &e) {
            if(rgba) {
                delete [] rgba;
                rgba = NULL;
            }
            if(colorBuffer) {
                delete [] colorBuffer;
                colorBuffer = NULL;
            }
            std::cerr << "unable to read iamge " << filename << 
                ":" << e.what();
            return NULL;
        }
    }

    static bool writeImageEXR(const string& filename, const Color* colorBuffer,
        int width, int height) {
        Rgba* hrgba = NULL;
        try {
            hrgba = new Rgba[width * height];
            for(int i = 0; i < width * height; ++i) {
                hrgba[i] = Rgba(
                    colorBuffer[i].r,
                    colorBuffer[i].g,
                    colorBuffer[i].b,
                    colorBuffer[i].a);
            }
            Header header(width, height);
            header.channels().insert("R", Channel(HALF));
            header.channels().insert("G", Channel(HALF));
            header.channels().insert("B", Channel(HALF));
            header.channels().insert("A", Channel(HALF));

            OutputFile file(filename.c_str(), header);

            FrameBuffer frameBuffer;
            frameBuffer.insert("R", Slice(
                HALF,                               // type
                (char*)(hrgba),                     // base
                4 * sizeof(half),                   // xstride 
                4 * sizeof(half) * width,           // ystride
                1, 1,                               // x/y sampling
                0.0));                              // fill value
            frameBuffer.insert("G", Slice(
                HALF,                               // type
                (char*)(hrgba) + 1 * sizeof(half),  // base
                4 * sizeof(half),                   // xstride 
                4 * sizeof(half) * width,           // ystride
                1, 1,                               // x/y sampling
                0.0));                              // fill value
            frameBuffer.insert("B", Slice(
                HALF,                               // type
                (char*)(hrgba) + 2 * sizeof(half),  // base
                4 * sizeof(half),                   // xstride 
                4 * sizeof(half) * width,           // ystride
                1, 1,                               // x/y sampling
                0.0));                              // fill value
            frameBuffer.insert("A", Slice(
                HALF,                               // type
                (char*)(hrgba) + 3 * sizeof(half),  // base
                4 * sizeof(half),                   // xstride 
                4 * sizeof(half) * width,           // ystride
                1, 1,                               // x/y sampling
                1.0));                              // fill value

            file.setFrameBuffer(frameBuffer);
            file.writePixels(height);
            delete [] hrgba;
            return true;
        } catch (const std::exception &e) {
            if(hrgba) {
                delete [] hrgba;
                hrgba = NULL;
            }
            std::cerr << "unable to read iamge " << filename << 
                ":" << e.what();
            return false;
        }
    }

#endif //GOBLIN_ENABLE_EXR

    static Color* loadImagePNG(const string& filename, int *width, int *height) {
        FILE* fp = fopen(filename.c_str(), "rb");
        if(!fp) {
            std::cerr << "can't open image file " << filename << std::endl;
            return NULL;
        } 

        // read in header
        const size_t HEADER_LENGTH = 8;
        png_byte header[HEADER_LENGTH];
        size_t n = fread(header, 1, HEADER_LENGTH, fp);
        if(n != HEADER_LENGTH || png_sig_cmp(header, 0, HEADER_LENGTH)) {
            std::cerr << "error: " << filename << " is not a png file\n";
            return NULL;
        }

        png_structp pngPtr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
            NULL, NULL, NULL);    
        if(!pngPtr) {
            fclose(fp);
            return NULL;
        }
        
        png_infop infoPtr = png_create_info_struct(pngPtr);
        if(!infoPtr) {
            png_destroy_read_struct(&pngPtr, (png_infopp)NULL, 
                (png_infopp)NULL);
            fclose(fp);
            return NULL;
        }

        png_infop endInfo = png_create_info_struct(pngPtr);
        if(!endInfo) {
            png_destroy_read_struct(&pngPtr, &infoPtr, (png_infopp)NULL);
            fclose(fp);
            return NULL;
        }
        
        if(setjmp(png_jmpbuf(pngPtr))) {
            png_destroy_read_struct(&pngPtr, &infoPtr, &endInfo);
            fclose(fp);
            return NULL;
        }

        png_init_io(pngPtr, fp);
        // let libpng know header bytes are already read
        png_set_sig_bytes(pngPtr, HEADER_LENGTH);
        
        // get the image info
        png_read_info(pngPtr, infoPtr);
        *width = png_get_image_width(pngPtr, infoPtr);
        *height = png_get_image_height(pngPtr, infoPtr);
        int bitDepth = png_get_bit_depth(pngPtr, infoPtr);
        png_byte colorType = png_get_color_type(pngPtr, infoPtr);

        // force the image info to RGBA 32bit
        if(colorType != PNG_COLOR_TYPE_RGBA) {
            png_set_expand(pngPtr);
        }
        if(colorType == PNG_COLOR_TYPE_GRAY ||
            colorType == PNG_COLOR_TYPE_GRAY_ALPHA) {
            png_set_gray_to_rgb(pngPtr);
        } 
        if(bitDepth < 8) {
            png_set_packing(pngPtr);
        } else if (bitDepth == 16) {
            png_set_strip_16(pngPtr);
        }
        if(colorType != PNG_COLOR_TYPE_RGBA) {
            png_set_filler(pngPtr, 255, PNG_FILLER_AFTER);
        }
        png_read_update_info(pngPtr, infoPtr);
        
        // make sure it's rgba 32 bit mode now
        if((int)png_get_rowbytes(pngPtr, infoPtr) != ((*width) * 4)) {
            png_destroy_read_struct(&pngPtr, &infoPtr, &endInfo);
            fclose(fp);
            return NULL;
        }
         
        // read in the file
        size_t bufferSize = (*width) * (*height) * 4;
        unsigned char* buffer = new unsigned char[bufferSize];
        png_bytep* rowPointers= new png_bytep[sizeof(png_bytep)* (*height)];
        for(int y = 0; y < (*height); ++y) {
            rowPointers[y] = &buffer[y * (*width) *4];
        }
        png_read_image(pngPtr, rowPointers);

        Color* colorBuffer = new Color[(*width) * (*height)];
        // fill in to new allocated Color buffer
        for(int y = 0; y < *height; ++y) {
            for(int x = 0; x < *width; ++x) {
                size_t index = y * (*width) + x;
                unsigned char* p = &buffer[index * 4];
                colorBuffer[index].r = (float)p[0] / 255.0f;
                colorBuffer[index].g = (float)p[1] / 255.0f;
                colorBuffer[index].b = (float)p[2] / 255.0f;
                colorBuffer[index].a = (float)p[3] / 255.0f;
            }
        }

        // clean up
        delete [] buffer;
        delete [] rowPointers; 
        fclose(fp);
        png_destroy_read_struct(&pngPtr, &infoPtr, &endInfo);
        return colorBuffer;
    }

   
    static bool writeImagePNG(const string& filename, const Color* colorBuffer,
            int width, int height, float gama = 2.2f) {
        FILE *fp = fopen(filename.c_str(), "wb");
        if(!fp) {
            return false;
        }

        png_structp pngPtr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                NULL, NULL, NULL);
        if(!pngPtr) {
            fclose(fp);
            return false;
        }
        png_infop infoPtr = png_create_info_struct(pngPtr);
        if(!infoPtr) {
            fclose(fp);
            png_destroy_write_struct(&pngPtr, (png_infopp)NULL);
            return false;
        }

        if(setjmp(png_jmpbuf(pngPtr))) {
            png_destroy_write_struct(&pngPtr, &infoPtr);
            fclose(fp);
            return false;
        }

        png_init_io(pngPtr, fp);

        // write out header
        png_set_IHDR(pngPtr, infoPtr, width , height, 8, 
            PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE, 
            PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
        png_write_info(pngPtr, infoPtr);

        //write out image data
        png_bytep *rowPointers = new png_bytep[height];
        unsigned char* buffer = new unsigned char[width * height * 4];
        for(int y = 0; y < height; ++y) {
            for(int x = 0; x < width; ++x) {
                size_t index = y * width + x;
                Color c = colorBuffer[index];
                // gama correction
                float invGama = 1.0f / gama;
                c.r = pow(c.r, invGama);
                c.g = pow(c.g, invGama);
                c.b = pow(c.b, invGama);
                buffer[index * 4 + 0] = 
                    static_cast<unsigned char>(clamp(c.r, 0, 1) * 255.0f);
                buffer[index * 4 + 1] =
                    static_cast<unsigned char>(clamp(c.g, 0, 1) * 255.0f);
                buffer[index * 4 + 2] = 
                    static_cast<unsigned char>(clamp(c.b, 0, 1) * 255.0f);
                buffer[index * 4 + 3] = 
                    static_cast<unsigned char>(clamp(c.a, 0, 1) * 255.0f);
            }
            rowPointers[y] = &buffer[y * width * 4];
        }
        png_write_image(pngPtr, rowPointers);
        png_write_end(pngPtr, infoPtr);
        // clean up
        delete [] rowPointers;
        delete [] buffer;
        fclose(fp);
        png_destroy_write_struct(&pngPtr, &infoPtr);
        return true;
    }


    Color* loadImage(const string& filename, int *width, int *height) {
        size_t extOffset = filename.rfind(".");

        if(extOffset == string::npos) {
            std::cerr << "error loading image " << filename <<
                " :unrecognized file format" << std::endl;
            return NULL;
        }
        string ext = filename.substr(extOffset);
        if(ext == ".png" || ext == ".PNG") {
            return loadImagePNG(filename, width, height);
#ifdef GOBLIN_ENABLE_EXR
        } else if(ext == ".exr" || ext == ".EXR") {
            return loadImageEXR(filename, width, height);
#endif //GOBLIN_ENABLE_EXR
        } else {
            // TODO .exr .tiff .tga
            std::cerr << "error loading image " << filename <<
                " :unsupported format " << std::endl;
            return NULL;
        }
    }

 
    bool writeImage(const string& filename, Color* colorBuffer,
            int width, int height, bool bloomImage) {

        if(bloomImage) {
            bloom(colorBuffer, width, height);
        }

        size_t extOffset = filename.rfind(".");
        if(extOffset == string::npos) {
            return writeImagePNG(filename + ".png", colorBuffer, width, height);
        }
        string ext = filename.substr(extOffset);
        if(ext == ".png" || ext == ".PNG") {
            toneMapping(colorBuffer, width, height);
            return writeImagePNG(filename, colorBuffer, width, height);
#ifdef GOBLIN_ENABLE_EXR
        } else if(ext == ".exr" || ext == ".EXR") {
            return writeImageEXR(filename, colorBuffer, width, height);
#endif //GOBLIN_ENABLE_EXR
        } else {
            // TODO .exr .tiff .tga
            std::cerr << "format " << ext << "is not supported yet" << 
                std::endl;
            return false;
        }
    }

    void bloom(Color* colorBuffer, int width, int height,
        float bloomRadius, float bloomWeight) {
        if(bloomRadius <= 0.0f || bloomWeight <= 0.0f) {
            return;
        }
        // build up filter lookup table
        int filterWidth = ceilInt(bloomRadius * max(width, height)) / 2;
        float* filter = new float[filterWidth * filterWidth];
        for(int y = 0; y < filterWidth; ++y) {
            for(int x = 0; x < filterWidth; ++x) {
                float d = sqrtf((float)(x * x + y * y)) / (float)filterWidth;
                filter[y * filterWidth + x] = powf(max(0.0f, 1.0f - d), 4.0f);
            }
        }

        Color* bloomResult = new Color[width * height];
        for(int i = 0; i < width * height; ++i) {
            *bloomResult = Color::Black;
        }
        // apply filter result to result layer
        for(int y = 0; y < height; ++y) {
            for(int x = 0; x < width; ++x) {
                int x0 = max(0, x - filterWidth + 1);
                int x1 = min(x + filterWidth - 1, width - 1);
                int y0 = max(0, y - filterWidth + 1);
                int y1 = min(y + filterWidth - 1, height - 1);
                int index = y * width + x;
                float weightSum = 0.0f;
                for(int py = y0; py <= y1; ++py) {
                    for(int px = x0; px <= x1; ++px) {
                        int fx = abs(px - x);
                        int fy = abs(py - y);
                        // only neighbor, exclude pixel(x, y) itself
                        if(fx == 0 && fy == 0) {
                            continue;
                        }
                        float weight = filter[fy * filterWidth + fx];
                        bloomResult[index] += 
                            weight * colorBuffer[py * width + px];
                        weightSum += weight;
                    }
                }
                bloomResult[index] /= weightSum;
            }
        }
        // blend input layer with bloomResult
        for(int i = 0; i < width * height; ++i) {
            colorBuffer[i] = (1.0f - bloomWeight) * colorBuffer[i] +
                bloomWeight * bloomResult[i];
        }
        delete [] filter;
        delete [] bloomResult;
    }

    // based on Reinhard. E Siggraph 02
    // Photographic Tone Reproduction for Digital Images
    void toneMapping(Color* colorBuffer, int width, int height) {
        // world adaptation luminance
        float Ywa = 0.0f;
        for(int i = 0; i < width * height; ++i) {
            float y = colorBuffer[i].luminance();
            Ywa += logf(1e4f + y);
        }
        Ywa = expf(Ywa / (width * height));
        float invy2 = 1.0f / (Ywa * Ywa);
        for(int i = 0; i < width * height; ++i) {
            float y = colorBuffer[i].luminance();
            float s = (1.0f + y * invy2) / (1.0f + y);
            colorBuffer[i] *= s;
        }
    }
}

