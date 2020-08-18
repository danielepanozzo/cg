// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "raster.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

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
		FragmentAttributes out(uniform.color(0),uniform.color(1),uniform.color(2),uniform.color(3));
		out.position = va.position;
		return out;
	};

	// The blending shader converts colors between 0 and 1 to uint8
	program.BlendingShader = [](const FragmentAttributes& fa, const FrameBufferAttributes& previous)
	{
		if (fa.position[2] < previous.depth)
		{
			FrameBufferAttributes out(fa.color[0]*255, fa.color[1]*255, fa.color[2]*255, fa.color[3]*255);
			out.depth = fa.position[2];
			return out;
		}
		else
			return previous;
	};

	// First triangle
	vector<VertexAttributes> vertices_1;
	vertices_1.push_back(VertexAttributes(-1,-1,0));
	vertices_1.push_back(VertexAttributes(1,-1,0));
	vertices_1.push_back(VertexAttributes(-1,1,0));

	// Second triangle
	vector<VertexAttributes> vertices_2;
	vertices_2.push_back(VertexAttributes(-1,-1,0.5));
	vertices_2.push_back(VertexAttributes(1,-1,0.5));
	vertices_2.push_back(VertexAttributes(1,1,0.5));

	// Draw first triangle in red
	uniform.color << 1,0,0,1;
	rasterize_triangles(program,uniform,vertices_1,frameBuffer);

	// Draw second triangle in blue
	uniform.color << 0,0,1,1;
	rasterize_triangles(program,uniform,vertices_2,frameBuffer);

	vector<uint8_t> image;
	framebuffer_to_uint8(frameBuffer,image);
	stbi_write_png("triangle.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows()*4);
	
	return 0;
}
