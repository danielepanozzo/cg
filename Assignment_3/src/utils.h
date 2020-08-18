#ifndef UTILS_H
#define UTILS_H

#include <Eigen/Dense>
#include <algorithm>
#include <vector>

unsigned char double_to_unsignedchar(const double d) {
	return round(std::max(std::min(1.,d),0.)*255);
}

void write_matrix_to_uint8(
	const Eigen::MatrixXd& R, const Eigen::MatrixXd& G,
	const Eigen::MatrixXd& B, const Eigen::MatrixXd& A,
	std::vector<uint8_t>& image)
{
	assert(R.rows() == G.rows() && G.rows() == B.rows() && B.rows() == A.rows());
	assert(R.cols() == G.cols() && G.cols() == B.cols() && B.cols() == A.cols());

	const int w = R.rows();                              // Image width
	const int h = R.cols();                              // Image height
	const int comp = 4;                                  // 4 Channels Red, Green, Blue, Alpha
	const int stride_in_bytes = w*comp;                  // Length of one row in bytes
	image.resize(w*h*comp,0);         // The image itself;

	for (unsigned wi = 0; wi < w; ++wi) {
		for (unsigned hi = 0; hi < h; ++hi) {
			image[(hi * w * 4) + (wi * 4) + 0] = double_to_unsignedchar(R(wi,hi));
			image[(hi * w * 4) + (wi * 4) + 1] = double_to_unsignedchar(G(wi,hi));
			image[(hi * w * 4) + (wi * 4) + 2] = double_to_unsignedchar(B(wi,hi));
			image[(hi * w * 4) + (wi * 4) + 3] = double_to_unsignedchar(A(wi,hi));
		}
	}
}

void write_matrix_to_png(
	const Eigen::MatrixXd& R, const Eigen::MatrixXd& G,
	const Eigen::MatrixXd& B, const Eigen::MatrixXd& A,
	const std::string& filename)
{
	const int w = R.rows();                              // Image width
	const int h = R.cols();                              // Image height
	const int comp = 4;                                  // 3 Channels Red, Green, Blue, Alpha
	const int stride_in_bytes = w*comp;                  // Length of one row in bytes

	std::vector<uint8_t> image;
	write_matrix_to_uint8(R,G,B,A,image);
	stbi_write_png(filename.c_str(), w, h, comp, image.data(), stride_in_bytes);

}

#endif
