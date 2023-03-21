#include "raster.h"	
#include <iostream>

void rasterize_triangle(const Program& program, const UniformAttributes& uniform, const VertexAttributes& v1, const VertexAttributes& v2, const VertexAttributes& v3, FrameBuffer& frameBuffer)
{
		// Collect coordinates into a matrix and convert to canonical representation
		Eigen::Matrix<double,3,4> p;
		p.row(0) = v1.position.array()/v1.position[3];
		p.row(1) = v2.position.array()/v2.position[3];
		p.row(2) = v3.position.array()/v3.position[3];

		// Coordinates are in -1..1, rescale to pixel size (x,y only)
		p.col(0) = ((p.col(0).array()+1.0)/2.0)*frameBuffer.rows();		
		p.col(1) = ((p.col(1).array()+1.0)/2.0)*frameBuffer.cols();

		// Find bounding box in pixels
		int lx = std::floor(p.col(0).minCoeff());
		int ly = std::floor(p.col(1).minCoeff());
		int ux = std::ceil(p.col(0).maxCoeff());
		int uy = std::ceil(p.col(1).maxCoeff());

		// Clamp to framebuffer
		lx = std::min(std::max(lx,int(0)),int(frameBuffer.rows()-1));
		ly = std::min(std::max(ly,int(0)),int(frameBuffer.cols()-1));
		ux = std::min(std::max(ux,int(0)),int(frameBuffer.rows()-1));
		uy = std::min(std::max(uy,int(0)),int(frameBuffer.cols()-1));

		// Build the implicit triangle representation
		Eigen::Matrix3d A;
		A.col(0) = p.row(0).segment(0,3);
		A.col(1) = p.row(1).segment(0,3);
		A.col(2) = p.row(2).segment(0,3);
		A.row(2) << 1.0, 1.0, 1.0;

		Eigen::Matrix3d Ai = A.inverse();

		// Rasterize the triangle
		for (unsigned i=lx; i<=ux; i++)
		{
			for (unsigned j=ly; j<=uy; j++)
			{
				// The pixel center is offset by 0.5, 0.5
				Eigen::Vector3d pixel(i+0.5,j+0.5,1);
				Eigen::Vector3d b = Ai*pixel;
				if (b.minCoeff() >= 0)
				{
					VertexAttributes va = VertexAttributes::interpolate(v1,v2,v3,b[0],b[1],b[2]);
					// Only render fragments within the bi-unit cube
					if (va.position[2] >= -1 && va.position[2] <= 1)
					{ 
						FragmentAttributes frag = program.FragmentShader(va,uniform);
						frameBuffer(i,j) = program.BlendingShader(frag,frameBuffer(i,j));
					}
				}
			}
		}
}

void rasterize_triangles(const Program& program, const UniformAttributes& uniform, const std::vector<VertexAttributes>& vertices, FrameBuffer& frameBuffer)
{
	// Call vertex shader on all vertices
	std::vector<VertexAttributes> v(vertices.size());
	for (unsigned i=0; i<vertices.size();i++)
		v[i] = program.VertexShader(vertices[i],uniform);

	// Call the rasterization function on every triangle
	for (unsigned i=0; i<vertices.size()/3; i++)
		rasterize_triangle(program,uniform,v[i*3+0],v[i*3+1],v[i*3+2],frameBuffer);
}

void rasterize_line(const Program& program, const UniformAttributes& uniform, const VertexAttributes& v1, const VertexAttributes& v2, double line_thickness, FrameBuffer& frameBuffer)
{
		// Collect coordinates into a matrix and convert to canonical representation
		Eigen::Matrix<double,2,4> p;
		p.row(0) = v1.position.array()/v1.position[3];
		p.row(1) = v2.position.array()/v2.position[3];

		// Coordinates are in -1..1, rescale to pixel size (x,y only)
		p.col(0) = ((p.col(0).array()+1.0)/2.0)*frameBuffer.rows();	
		p.col(1) = ((p.col(1).array()+1.0)/2.0)*frameBuffer.cols();

		// Find bounding box in pixels, adding the line thickness
		int lx = std::floor(p.col(0).minCoeff()-line_thickness);
		int ly = std::floor(p.col(1).minCoeff()-line_thickness);
		int ux = std::ceil(p.col(0).maxCoeff()+line_thickness);
		int uy = std::ceil(p.col(1).maxCoeff()+line_thickness);

		// Clamp to framebuffer
		lx = std::min(std::max(lx,int(0)),int(frameBuffer.rows()-1));
		ly = std::min(std::max(ly,int(0)),int(frameBuffer.cols()-1));
		ux = std::min(std::max(ux,int(0)),int(frameBuffer.rows()-1));
		uy = std::min(std::max(uy,int(0)),int(frameBuffer.cols()-1));

		// We only need the 2d coordinates of the endpoints of the line
		Eigen::Vector2f l1(p(0,0),p(0,1));
		Eigen::Vector2f l2(p(1,0),p(1,1));

		// Parametrize the line as l1 + t (l2-l1)
		double t = -1;
		double ll  = (l1-l2).squaredNorm();

		// Rasterize the line
		for (unsigned i=lx; i<=ux; i++)
		{
			for (unsigned j=ly; j<=uy; j++)
			{
				// The pixel center is offset by 0.5, 0.5
				Eigen::Vector2f pixel(i+0.5,j+0.5);

				if (ll == 0.0)
					// The segment has zero length
					t = 0;
				else
				{
					// Project p on the line
					t = (pixel-l1).dot(l2-l1)/ll;
					// Clamp between 0 and 1
					t = std::fmax(0, std::fmin(1, t));
				}

  				Eigen::Vector2f pixel_p = l1 + t * (l2 - l1);
				
				if ((pixel - pixel_p).squaredNorm() < (line_thickness*line_thickness))
				{
					VertexAttributes va = VertexAttributes::interpolate(v1,v2,v1,1-t,t,0);
					FragmentAttributes frag = program.FragmentShader(va,uniform);
					frameBuffer(i,j) = program.BlendingShader(frag,frameBuffer(i,j));
				}
			}
		}
}

void rasterize_lines(const Program& program, const UniformAttributes& uniform, const std::vector<VertexAttributes>& vertices, double line_thickness, FrameBuffer& frameBuffer)
{
	// Call vertex shader on all vertices
	std::vector<VertexAttributes> v(vertices.size());
	for (unsigned i=0; i<vertices.size();i++)
		v[i] = program.VertexShader(vertices[i],uniform);

	// Call the rasterization function on every line
	for (unsigned i=0; i<vertices.size()/2; i++)
		rasterize_line(program,uniform,v[i*2+0],v[i*2+1],line_thickness,frameBuffer);
}

void framebuffer_to_uint8(const FrameBuffer& frameBuffer, std::vector<uint8_t>& image)
{
	const int w = frameBuffer.rows();                              // Image width
	const int h = frameBuffer.cols();                              // Image height
	const int comp = 4;                                  // 4 Channels Red, Green, Blue, Alpha
	const int stride_in_bytes = w*comp;                  // Length of one row in bytes
	image.resize(w*h*comp,0);         // The image itself;

	for (unsigned wi = 0; wi < w; ++wi) 
	{
		for (unsigned hi = 0; hi < h; ++hi) 
		{
			unsigned hif = h-1-hi;
			image[(hi * w * 4) + (wi * 4) + 0] = frameBuffer(wi,hif).color[0];
			image[(hi * w * 4) + (wi * 4) + 1] = frameBuffer(wi,hif).color[1];
			image[(hi * w * 4) + (wi * 4) + 2] = frameBuffer(wi,hif).color[2];
			image[(hi * w * 4) + (wi * 4) + 3] = frameBuffer(wi,hif).color[3];
		}
	}
}
