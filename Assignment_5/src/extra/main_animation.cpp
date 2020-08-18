// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "raster.h"
#include <gif.h>

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include <stb_image_write.h>

using namespace std;

int main() 
{

	// The Framebuffer storing the image rendered by the rasterizer
	Eigen::Matrix<FrameBufferAttributes,Eigen::Dynamic,Eigen::Dynamic> frameBuffer(500,500);

	// Global Constants (empty in this example)
	UniformAttributes uniform;

	// Basic rasterization program
	Program program;

	// The vertex shader is the identity
	program.VertexShader = [](const VertexAttributes& va, const UniformAttributes& uniform)
	{
		return va;
	};

	// The fragment shader uses a fixed color
	program.FragmentShader = [](const VertexAttributes& va, const UniformAttributes& uniform)
	{
		return FragmentAttributes(1,0,0);
	};

	// The blending shader converts colors between 0 and 1 to uint8
	program.BlendingShader = [](const FragmentAttributes& fa, const FrameBufferAttributes& previous)
	{
		return FrameBufferAttributes(fa.color[0]*255,fa.color[1]*255,fa.color[2]*255,fa.color[3]*255);
	};

	// One triangle in the center of the screen
	vector<VertexAttributes> vertices;
	vertices.push_back(VertexAttributes(-1,-1,0));
	vertices.push_back(VertexAttributes(1,-1,0));
	vertices.push_back(VertexAttributes(0,1,0));

	const char * fileName = "triangle.gif";
	vector<uint8_t> image;
	int delay = 25;
	GifWriter g;
	GifBegin(&g, fileName, frameBuffer.rows(), frameBuffer.cols(), delay);

	for (float i=0;i<1;i+=0.05)
	{
		frameBuffer.setConstant(FrameBufferAttributes());
		vertices[2].position[1] -= 0.05;
		rasterize_triangles(program,uniform,vertices,frameBuffer);
		framebuffer_to_uint8(frameBuffer,image);
		GifWriteFrame(&g, image.data(), frameBuffer.rows(), frameBuffer.cols(), delay);
	}

	GifEnd(&g);
	return 0;
}
