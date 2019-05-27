#include "GoblinImageIO.h"
#include "GoblinColor.h"
#include "GoblinUtils.h"

#include <iostream>
#include <cstdio>
#include <cassert>

#define TINYEXR_IMPLEMENTATION
#include "tinyexr.h"

namespace Goblin {

static Color* loadImageEXR(const string& filename,
    int *width, int *height) {
	static_assert(sizeof(Color) == 4 * sizeof(float),
		"loadImageEXR assumes class Color in float4 rgba memory format");
	float* out = nullptr;
	const char* err = nullptr;
	int ret = LoadEXR(&out, width, height, filename.c_str(), &err);
	if (ret != TINYEXR_SUCCESS) {
		std::cerr << "unable to read iamge " << filename;
		if (err) {
			std::cerr << " : " << err << std::endl;
			FreeEXRErrorMessage(err); // release memory of error message.
		} else {
			std::cerr << std::endl;
		}
		return nullptr;
	} else {
		return reinterpret_cast<Color*>(out);
	}
}

static bool writeImageEXR(const string& filename, const Color* colorBuffer,
    int width, int height) {

	EXRHeader header;
	InitEXRHeader(&header);

	EXRImage image;
	InitEXRImage(&image);

	image.num_channels = 3;

	std::vector<float> images[3];
	size_t numPixels = (size_t)width * (size_t)height;
	images[0].resize(numPixels);
	images[1].resize(numPixels);
	images[2].resize(numPixels);

	for (int i = 0; i < width * height; i++) {
		images[0][i] = colorBuffer[i].r;
		images[1][i] = colorBuffer[i].g;
		images[2][i] = colorBuffer[i].b;
	}

	float* imagePtr[3];
	imagePtr[0] = &(images[2].at(0)); // B
	imagePtr[1] = &(images[1].at(0)); // G
	imagePtr[2] = &(images[0].at(0)); // R

	image.images = (unsigned char**)imagePtr;
	image.width = width;
	image.height = height;

	header.num_channels = 3;
	header.channels = (EXRChannelInfo*)malloc(
		sizeof(EXRChannelInfo) * header.num_channels);
	// Must be BGR(A) order, since most of EXR viewers expect this channel order.
	header.channels[0].name[0] = 'B';
	header.channels[0].name[1] = '\0';
	header.channels[1].name[0] = 'G';
	header.channels[1].name[1] = '\0';
	header.channels[2].name[0] = 'R';
	header.channels[2].name[1] = '\0';

	header.pixel_types = (int*)malloc(sizeof(int) * header.num_channels);
	header.requested_pixel_types = (int*)malloc(sizeof(int) * header.num_channels);
	for (int i = 0; i < header.num_channels; i++) {
		// pixel type of input image
		header.pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT;
		// pixel type of output image to be stored in .EXR
		header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_HALF;
	}

	const char* err;
	int ret = SaveEXRImageToFile(&image, &header, filename.c_str(), &err);
	free(header.channels);
	free(header.pixel_types);
	free(header.requested_pixel_types);
	if (ret != TINYEXR_SUCCESS) {
		std::cerr << "unable to write iamge " << filename << " : " <<
			err << std::endl;
		return false;
	}
	return true;
}

static bool writeImagePPM(const string& filename, const Color* colorBuffer,
	int width, int height, float gamma = 2.2f) {
	FILE* fp;
	errno_t err = fopen_s(&fp, filename.c_str(), "w");
	if (err != 0 || fp == nullptr) {
		std::cerr << "can not open file " << filename << std::endl;
		return false;
	}

	fprintf(fp, "P3\n%d %d\n%d\n", width, height, 255);
	float invGama = 1.0f / gamma;
	for (int y = 0; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
			int index = y * width + x;
			Color c = colorBuffer[index];
			// gama correction
			c.r = pow(c.r, invGama);
			c.g = pow(c.g, invGama);
			c.b = pow(c.b, invGama);
			fprintf(fp, "%d %d %d ",
				static_cast<int>(clamp(c.r, 0.0f, 1.0f) * 255.0f),
				static_cast<int>(clamp(c.g, 0.0f, 1.0f) * 255.0f),
				static_cast<int>(clamp(c.b, 0.0f, 1.0f) * 255.0f));
		}
	}
	return true;
}

Color* loadImage(const string& filename, int *width, int *height) {
    size_t extOffset = filename.rfind(".");

    if (extOffset == string::npos) {
        std::cerr << "error loading image " << filename <<
            " :unrecognized file format" << std::endl;
        return NULL;
    }
    string ext = filename.substr(extOffset);
    if (ext == ".exr" || ext == ".EXR") {
        return loadImageEXR(filename, width, height);
    } else {
        // TODO .tiff .tga
        std::cerr << "error loading image " << filename <<
            " :unsupported format " << ext << std::endl;
        return NULL;
    }
}
 
bool writeImage(const string& filename, Color* colorBuffer,
        int width, int height, bool doToneMapping) {
    size_t extOffset = filename.rfind(".");
    if (extOffset == string::npos) {
        return writeImagePPM(filename + ".ppm", colorBuffer, width, height);
    }
    string ext = filename.substr(extOffset);
    if (ext == ".ppm" || ext == ".PPM") {
        if (doToneMapping) {
            toneMapping(colorBuffer, width, height);
        }
        return writeImagePPM(filename, colorBuffer, width, height);
    } else if (ext == ".exr" || ext == ".EXR") {
        return writeImageEXR(filename, colorBuffer, width, height);
    } else {
        // TODO .tiff .tga
        std::cerr << "format " << ext << "is not supported yet" <<
            std::endl;
        return writeImagePPM(filename + ".ppm", colorBuffer, width, height);
    }
}

void bloom(Color* colorBuffer, int width, int height,
    float bloomRadius, float bloomWeight) {
    if (bloomRadius <= 0.0f || bloomWeight <= 0.0f) {
        return;
    }
    // build up filter lookup table
    int filterWidth = ceilInt(bloomRadius * max(width, height)) / 2;
    float* filter = new float[filterWidth * filterWidth];
    for (int y = 0; y < filterWidth; ++y) {
        for (int x = 0; x < filterWidth; ++x) {
            float d = sqrtf((float)(x * x + y * y)) / (float)filterWidth;
            filter[y * filterWidth + x] = powf(max(0.0f, 1.0f - d), 4.0f);
        }
    }

    Color* bloomResult = new Color[width * height];
    for (int i = 0; i < width * height; ++i) {
        bloomResult[i] = Color::Black;
    }
    // apply filter result to result layer
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            int x0 = max(0, x - filterWidth + 1);
            int x1 = min(x + filterWidth - 1, width - 1);
            int y0 = max(0, y - filterWidth + 1);
            int y1 = min(y + filterWidth - 1, height - 1);
            int index = y * width + x;
            float weightSum = 0.0f;
            for (int py = y0; py <= y1; ++py) {
                for (int px = x0; px <= x1; ++px) {
                    int fx = abs(px - x);
                    int fy = abs(py - y);
                    // only neighbor, exclude pixel(x, y) itself
                    if (fx == 0 && fy == 0) {
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
    for (int i = 0; i < width * height; ++i) {
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
    for (int i = 0; i < width * height; ++i) {
        float y = colorBuffer[i].luminance();
        Ywa += logf(1e4f + y);
    }
    Ywa = expf(Ywa / (width * height));
    float invy2 = 1.0f / (Ywa * Ywa);
    for (int i = 0; i < width * height; ++i) {
        float y = colorBuffer[i].luminance();
        float s = (1.0f + y * invy2) / (1.0f + y);
        colorBuffer[i] *= s;
    }
}

}