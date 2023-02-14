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

			// Generic implementation of ray sphere intersection
            // Based on class notes, formulating intersection as (d.d)t^2 + 2d.(e-c)t + (e-c)(e-c) - R^2 = 0 where e is origin of ray and d is direction of ray
            double _A = ray_direction.dot(ray_direction);
            double _B = 2 * ray_direction.dot(ray_origin - sphere_origin);
            double _C = (ray_origin-sphere_origin).dot(ray_origin-sphere_origin) - sphere_radius*sphere_radius;

            double discriminant = _B*_B - 4*_A*_C;

            if (discriminant >= 0) {
                double t1 = (-_B + pow(discriminant, 0.5)) / (2*_A);
                double t2 = (-_B - pow(discriminant, 0.5)) / (2*_A);

                double t_int = (t1 > t2) ? t2: t1;
                Vector3d ray_intersection = ray_origin + t_int * ray_direction;
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
	Vector3d pgram_origin(-0.75, -0.75, 0);
	Vector3d pgram_u(1, 0, 1);
	Vector3d pgram_v(0.5, 0.866, 0);

	// Single light source
	const Vector3d light_position(-1,1,1);
    Vector3d camera_position(0, 2, 5);

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

void raytrace_shading(){
	std::cout << "Simple ray tracer, one sphere with different shading" << std::endl;

	const std::string filename("shading.png");
	MatrixXd C = MatrixXd::Zero(800,800); // Store the color
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1,1,1);
	Vector3d x_displacement(2.0/C.cols(),0,0);
	Vector3d y_displacement(0,-2.0/C.rows(),0);

	// Single light source
	const Vector3d light_position(-1,1,1);
	double ambient = 0.1;
	MatrixXd diffuse = MatrixXd::Zero(800, 800);
	MatrixXd specular = MatrixXd::Zero(800, 800);

	for (unsigned i=0; i < C.cols(); ++i) {
		for (unsigned j=0; j < C.rows(); ++j) {
			// Prepare the ray
			Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d ray_direction = RowVector3d(0,0,-1);

			// Intersect with the sphere
			// NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
			Vector2d ray_on_xy(ray_origin(0),ray_origin(1));
			const double sphere_radius = 0.9;

			if (ray_on_xy.norm() < sphere_radius) {
				// The ray hit the sphere, compute the exact intersection point
				Vector3d ray_intersection(ray_on_xy(0),ray_on_xy(1),sqrt(sphere_radius*sphere_radius - ray_on_xy.squaredNorm()));

				// Compute normal at the intersection point
				Vector3d ray_normal = ray_intersection.normalized();

				// TODO: Add shading parameter here
				diffuse(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;
				specular(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;

				// Simple diffuse model
				C(i,j) = ambient + diffuse(i,j) + specular(i,j);

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

int main() {
//	raytrace_sphere();
//	raytrace_parallelogram();
	raytrace_perspective();
//	raytrace_shading();

	return 0;
}
