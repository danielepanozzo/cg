// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"


// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

bool check_ray_sphere_intersection(Vector3d& ray_origin, Vector3d& ray_direction, Vector3d& sphere_origin,
                                   double sphere_radius, double& t) {
    // Based on class notes, formulating intersection as (d.d)t^2 + 2d.(e-c)t + (e-c)(e-c) - R^2 = 0 where e is origin of ray and d is direction of ray
    double _A = ray_direction.dot(ray_direction);
    double _B = 2 * ray_direction.dot(ray_origin - sphere_origin);
    double _C = (ray_origin-sphere_origin).dot(ray_origin-sphere_origin) - sphere_radius*sphere_radius;

    double discriminant = _B*_B - 4*_A*_C;
    if (discriminant < 0) {
        return false;
    }

    double t1 = (-_B + pow(discriminant, 0.5)) / (2*_A);
    double t2 = (-_B - pow(discriminant, 0.5)) / (2*_A);
    t = (t1 > t2) ? t2: t1;
    return true;
}

void raytrace_sphere() {
	std::cout << "Simple ray tracer, one sphere with orthographic projection" << std::endl;

	const std::string filename("sphere_orthographic.png");
	MatrixXd C = MatrixXd::Zero(800,800); // Store the color
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1,1,1);
	Vector3d x_displacement(2.0/C.cols(),0,0);
	Vector3d y_displacement(0,-2.0/C.rows(),0);

	// Single light source
	const Vector3d light_position(-1,1,1);

    Vector3d sphere_origin(0, 0, 0);
    const double sphere_radius = 0.9;
    Vector3d w(0, 0, 1); // orthographic projection image plane normal direction

	for (unsigned i=0; i < C.cols(); ++i) {
		for (unsigned j=0; j < C.rows(); ++j) {
			// Prepare the ray
			Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d ray_direction = -1 * w;

            double t = 0;
            if (check_ray_sphere_intersection(ray_origin, ray_direction, sphere_origin, sphere_radius, t)) {
                Vector3d ray_intersection = ray_origin + t * ray_direction;
                Vector3d ray_normal = (ray_intersection-sphere_origin).normalized();
                C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;
                C(i,j) = std::max(C(i,j),0.);
                A(i,j) = 1;
            }
		}
	}

	// Save to png
	write_matrix_to_png(C,C,C,A,filename);

}

bool check_ray_parallelogram_intersection(Vector3d& ray_origin, Vector3d& ray_direction, Vector3d& pgram_origin,
                                          Vector3d& pgram_u, Vector3d& pgram_v, double& t) {
    Matrix3d A;
    Vector3d b;
    A << pgram_u.x(), pgram_v.x(), -ray_direction.x(), pgram_u.y(), pgram_v.y(), -ray_direction.y(), pgram_u.z(), pgram_v.z(), -ray_direction.z();
    b = ray_origin - pgram_origin;

    Vector3d x = A.colPivHouseholderQr().solve(b);
    double u = x.x();
    double v = x.y();
    t = x.z();
    return u>=0 && u<=1 && v>=0 && v<=1 && t>0;
}

void raytrace_parallelogram() {
	std::cout << "Simple ray tracer, one parallelogram with orthographic projection" << std::endl;

	const std::string filename("plane_orthographic.png");
	MatrixXd C = MatrixXd::Zero(800,800); // Store the color
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1,1,1);
	Vector3d x_displacement(2.0/C.cols(),0,0);
	Vector3d y_displacement(0,-2.0/C.rows(),0);

	// Parameters of the parallelogram (position of the lower-left corner + two sides)
	Vector3d pgram_origin(-0.75, -0.75, 0);
	Vector3d pgram_u(1, 0, 0);
	Vector3d pgram_v(0.5, 0.866, 0);

	// Single light source
	const Vector3d light_position(-1,1,1);

	for (unsigned i=0; i < C.cols(); ++i) {
		for (unsigned j=0; j < C.rows(); ++j) {
			// Prepare the ray
			Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d ray_direction = RowVector3d(0,0,-1);

            double t = 0;
			if (check_ray_parallelogram_intersection(ray_origin, ray_direction, pgram_origin, pgram_u, pgram_v, t)) {
				// The ray hit the parallelogram, compute the exact intersection point
				Vector3d ray_intersection = ray_origin + t * ray_direction;

				// Compute normal at the intersection point
				Vector3d ray_normal = pgram_u.cross(pgram_v).normalized();
                if (ray_normal.dot(light_position-ray_intersection) < 0) {
                    ray_normal = -1 * ray_normal;
                }

				// Simple diffuse model
				C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;

				// Clamp to zero
				C(i,j) = std::max(C(i,j),0.);

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}
		}
	}

	// Save to png
	write_matrix_to_png(C,C,C,A,filename);
}

void raytrace_perspective() {
	std::cout << "Simple ray tracer, one parallelogram with perspective projection" << std::endl;

	const std::string filename("plane_perspective.png");
	MatrixXd C = MatrixXd::Zero(800,800); // Store the color
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1,1,1);
	Vector3d x_displacement(2.0/C.cols(),0,0);
	Vector3d y_displacement(0,-2.0/C.rows(),0);

	// Parameters of the parallelogram (position of the lower-left corner + two sides)
	Vector3d pgram_origin(-0.75, 0.4, 0);
	Vector3d pgram_u(0.7, 0.7, 0);
	Vector3d pgram_v(1, 0, 0);

	// Single light source
	const Vector3d light_position(-1,1,1);
    Vector3d camera_position(0, 0, 5);
    Vector3d sphere_origin(0, 0, 0);
    const double sphere_radius = 0.3;

	for (unsigned i=0; i < C.cols(); ++i) {
		for (unsigned j=0; j < C.rows(); ++j) {
			// Prepare the ray (origin point and direction)
			Vector3d ray_origin = camera_position;
			Vector3d ray_direction = (origin + double(i)*x_displacement + double(j)*y_displacement) - camera_position;

            double t;
			// Check if the ray intersects with the parallelogram
			if (check_ray_parallelogram_intersection(ray_origin, ray_direction, pgram_origin, pgram_u, pgram_v, t)) {
				// The ray hit the parallelogram, compute the exact intersection point
				Vector3d ray_intersection = ray_origin + t * ray_direction;

				// Compute normal at the intersection point
				Vector3d ray_normal = pgram_u.cross(pgram_v).normalized();
                if (ray_normal.dot(light_position-ray_intersection) < 0) {
                    ray_normal = -1 * ray_normal;
                }

				// Simple diffuse model
				C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;

				// Clamp to zero
				C(i,j) = std::max(C(i,j),0.);

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}

            if (check_ray_sphere_intersection(ray_origin, ray_direction, sphere_origin, sphere_radius, t)) {
                Vector3d ray_intersection = ray_origin + t * ray_direction;
                Vector3d ray_normal = (ray_intersection-sphere_origin).normalized();
                C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;
                C(i,j) = std::max(C(i,j),0.);
                A(i,j) = 1;
            }
		}
	}

	// Save to png
	write_matrix_to_png(C,C,C,A,filename);
}

void raytrace_shading(){
	std::cout << "Simple ray tracer, one sphere with different shading" << std::endl;

	const std::string filename("shading.png");
	MatrixXd R = MatrixXd::Zero(800,800);
    MatrixXd G = MatrixXd::Zero(800,800);
    MatrixXd B = MatrixXd::Zero(800,800);
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1,1,1);
	Vector3d x_displacement(2.0/R.cols(),0,0);
	Vector3d y_displacement(0,-2.0/R.rows(),0);

	// Single light source
	const Vector3d light_position(1,1,1);
	double ambientR = 0.1;
    double ambientG = 0.1;
    double ambientB = 0.3;
    Vector3d diffuse_coeff(0.1, 0.1, 0.7);
    Vector3d specular_coeff(0.1, 0.1, 0.5);
    double phong_exp = 150;
	MatrixXd diffuseR = MatrixXd::Zero(800, 800);
    MatrixXd diffuseG = MatrixXd::Zero(800, 800);
    MatrixXd diffuseB = MatrixXd::Zero(800, 800);
    MatrixXd specularR = MatrixXd::Zero(800, 800);
    MatrixXd specularG = MatrixXd::Zero(800, 800);
	MatrixXd specularB = MatrixXd::Zero(800, 800);

    Vector3d sphere_origin(0, 0, 0);
    const double sphere_radius = 0.9;
    Vector3d w(0, 0, 1); // orthographic projection image plane normal direction

	for (unsigned i=0; i < R.cols(); ++i) {
		for (unsigned j=0; j < R.rows(); ++j) {
			// Prepare the ray
			Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d ray_direction = -1 * w;

            double t;
            if (check_ray_sphere_intersection(ray_origin, ray_direction, sphere_origin, sphere_radius, t)) {
                Vector3d ray_intersection = ray_origin + t * ray_direction;
                Vector3d ray_normal = (ray_intersection-sphere_origin).normalized();

                double diffuse = (light_position-ray_intersection).normalized().transpose() * ray_normal;
                diffuseR(i,j) = diffuse_coeff(0) * diffuse;
                diffuseG(i, j) = diffuse_coeff(1) * diffuse;
                diffuseB(i, j) = diffuse_coeff(2) * diffuse;

                Vector3d h = ((light_position-ray_intersection).normalized() + (ray_origin-ray_intersection).normalized()).normalized();
                double specular = pow(std::max(0., ray_normal.dot(h)), phong_exp);
                specularR(i,j) = specular_coeff(0) * specular;
                specularG(i,j) = specular_coeff(1) * specular;
                specularB(i,j) = specular_coeff(2) * specular;

                R(i,j) = ambientR + diffuseR(i,j) + specularR(i,j);
                R(i,j) = std::max(0., R(i,j));
                G(i,j) = ambientG + diffuseG(i,j) + specularG(i,j);
                G(i,j) = std::max(0., G(i,j));
                B(i,j) = ambientB + diffuseB(i,j) + specularB(i,j);
                B(i,j) = std::max(0., B(i,j));
                A(i,j) = 1;
            }
		}
	}

	// Save to png
	write_matrix_to_png(R,G,B,A,filename);
}

int main() {
//	raytrace_sphere();
//	raytrace_parallelogram();
	raytrace_perspective();
//	raytrace_shading();

	return 0;
}
