#include "GoblinImageIO.h"
#include "GoblinColor.h"
#include "GoblinUtils.h"

#include <iostream>
#include <cstdio>
#include <cassert>
#include <png.h>

namespace Goblin {
    static Color* loadImagePNG(const std::string& filename, int *width, int *height) {
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

   
    static bool writeImagePNG(const std::string& filename, const Color* colorBuffer,
            int width, int height) {
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

        // TODO 
        // add in tone mapping / gama correction to make this float->byte
        // conversion......less wasteful

        //write out image data
        png_bytep *rowPointers = new png_bytep[height];
        unsigned char* buffer = new unsigned char[width * height * 4];
        for(int y = 0; y < height; ++y) {
            for(int x = 0; x < width; ++x) {
                size_t index = y * width + x;
                buffer[index * 4 + 0] = 
                    static_cast<unsigned char>(clamp(colorBuffer[index].r, 0, 1) * 255.0f);
                buffer[index * 4 + 1] =
                    static_cast<unsigned char>(clamp(colorBuffer[index].g, 0, 1) * 255.0f);
                buffer[index * 4 + 2] = 
                    static_cast<unsigned char>(clamp(colorBuffer[index].b, 0, 1) * 255.0f);
                buffer[index * 4 + 3] = 
                    static_cast<unsigned char>(clamp(colorBuffer[index].a, 0, 1) * 255.0f);
            }
            rowPointers[y] = &buffer[y * width * 4];
        }
        png_write_image(pngPtr, rowPointers);
        png_write_end(pngPtr, infoPtr);

        //TODO replace this constant new/delete with better 
        //memory management pool
        
        // clean up
        delete [] rowPointers;
        delete [] buffer;
        fclose(fp);
        png_destroy_write_struct(&pngPtr, &infoPtr);
        return true;
    }


    Color* loadImage(const std::string& filename, int *width, int *height) {
        size_t extOffset = filename.rfind(".");

        if(extOffset == std::string::npos) {
            std::cerr << "error loading image " << filename <<
                " :unrecognized file format" << std::endl;
            return NULL;
        }
        std::string ext = filename.substr(extOffset);
        if(ext == ".png" || ext == ".PNG") {
            return loadImagePNG(filename, width, height);
        } else {
            // TODO .exr .tiff .tga
            std::cerr << "error loading image " << filename <<
                " :unsupported format " << std::endl;
            return NULL;
        }
    }

 
    bool writeImage(const std::string& filename, const Color* colorBuffer,
            int width, int height) {

        size_t extOffset = filename.rfind(".");
        if(extOffset == std::string::npos) {
            return writeImagePNG(filename + ".png", colorBuffer, width, height);
        }
        std::string ext = filename.substr(extOffset);
        if(ext == ".png" || ext == ".PNG") {
            return writeImagePNG(filename, colorBuffer, width, height);
        } else {
            // TODO .exr .tiff .tga
            std::cerr << "format " << ext << "is not supported yet" << 
                std::endl;
            return false;
        }
    }

}

