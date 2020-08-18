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
		VertexAttributes out;
		out.position = uniform.projective * va.position;
		return out;
	};

	// The fragment shader uses a fixed color
	program.FragmentShader = [](const VertexAttributes& va, const UniformAttributes& uniform)
	{
		return FragmentAttributes(1,1,1);
	};

	// The blending shader converts colors between 0 and 1 to uint8
	program.BlendingShader = [](const FragmentAttributes& fa, const FrameBufferAttributes& previous)
	{
		return FrameBufferAttributes(fa.color[0]*255,fa.color[1]*255,fa.color[2]*255,fa.color[3]*255);
	};

	// Builds a wireframe cube
	vector<VertexAttributes> vc;

	vc.push_back(VertexAttributes(0.75, 0.75, 0.25));
	vc.push_back(VertexAttributes(0.75, -0.75, 0.25));
	vc.push_back(VertexAttributes(0.75, 0.75, 0.75));
	vc.push_back(VertexAttributes(0.75, -0.75, 0.75));
	vc.push_back(VertexAttributes(-0.75, 0.75, 0.25));
	vc.push_back(VertexAttributes(-0.75, -0.75, 0.25));
	vc.push_back(VertexAttributes(-0.75, 0.75, 0.75));
	vc.push_back(VertexAttributes(-0.75, -0.75, 0.75));

	vector<VertexAttributes> vertices;

	vertices.push_back(vc[5]);
	vertices.push_back(vc[7]);

	vertices.push_back(vc[1]);
	vertices.push_back(vc[5]);

	vertices.push_back(vc[0]);
	vertices.push_back(vc[1]);

	vertices.push_back(vc[7]);
	vertices.push_back(vc[6]);

	vertices.push_back(vc[2]);
	vertices.push_back(vc[3]);

	vertices.push_back(vc[4]);
	vertices.push_back(vc[5]);

	vertices.push_back(vc[2]);
	vertices.push_back(vc[6]);

	vertices.push_back(vc[0]);
	vertices.push_back(vc[2]);

	vertices.push_back(vc[7]);
	vertices.push_back(vc[3]);

	vertices.push_back(vc[6]);
	vertices.push_back(vc[4]);

	vertices.push_back(vc[4]);
	vertices.push_back(vc[0]);

	vertices.push_back(vc[3]);
	vertices.push_back(vc[1]);

	// Add a projective transformation
	uniform.projective <<
	1, 0, 0, 0,
	0, 1, 0, 0,
	0, 0, 1, 0,
	0, 0, 1, 1;	

	rasterize_lines(program,uniform,vertices,2,frameBuffer);

	vector<uint8_t> image;
	framebuffer_to_uint8(frameBuffer,image);
	stbi_write_png("triangle.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows()*4);
	
	return 0;
}
