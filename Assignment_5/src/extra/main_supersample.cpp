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

FrameBuffer scale_down_4x(const FrameBuffer& fb)
{
	// The size of the framebuffer must be a multiple of 4
	assert(fb.rows() % 4 == 0);
	assert(fb.cols() % 4 == 0);

	// Allocate the reduced buffer
	FrameBuffer out(fb.rows()/4,fb.cols()/4);

	for (unsigned i=0;i<out.rows();i++)
	{
		for (unsigned j=0;j<out.cols();j++)
		{
			Eigen::Vector4f avg = Eigen::Vector4f::Zero();
			for (unsigned ii=0;ii<4;ii++)
				for (unsigned jj=0;jj<4;jj++)
					avg += fb(i*4+ii,j*4+jj).color.cast<float>();
			avg /= 16;
			out(i,j).color = avg.cast<uint8_t>();
		}
	}
	return out;
}

int main() 
{

	// The Framebuffer storing the image rendered by the rasterizer
	Eigen::Matrix<FrameBufferAttributes,Eigen::Dynamic,Eigen::Dynamic> frameBuffer(40*4,40*4);

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

	rasterize_triangles(program,uniform,vertices,frameBuffer);

	vector<uint8_t> image;
	frameBuffer = scale_down_4x(frameBuffer);
	framebuffer_to_uint8(frameBuffer,image);
	stbi_write_png("triangle.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows()*4);
	
	return 0;
}
