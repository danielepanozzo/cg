////////////////////////////////////////////////////////////////////////////////
// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <random>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

////////////////////////////////////////////////////////////////////////////////
// Scene setup, global variables
////////////////////////////////////////////////////////////////////////////////

const std::string filename("raytrace.png");

//Camera settings
// DOF changes
// Focal length set to focus on sphere 3 (just right of the checker colored sphere),
// z=-1 for this sphere and camera is at z=5, therefore focal_length=5-(-1)=6 should focus on this.
// Adjusted image_z to 1 (focal_length - camera's z) since image plane is at -image_z
const double focal_length = 6;
const double field_of_view = 0.7854; //45 degrees
const double image_z = 1;
const bool is_perspective = true;
const Vector3d camera_position(0, 0, 5);

const double epsilon = 0.00001; // shoot ray from epsilon away from source towards light
//Maximum number of recursive calls
const int max_bounce = 5;

// Objects
std::vector<Vector3d> sphere_centers;
std::vector<double> sphere_radii;
std::vector<Matrix3d> parallelograms;

//Material for the object, same material for all objects
const Vector4d obj_ambient_color(0.5, 0.1, 0.1, 0);
const Vector4d obj_diffuse_color(0.5, 0.5, 0.5, 0);
const Vector4d obj_specular_color(0.2, 0.2, 0.2, 0);
const double obj_specular_exponent = 256.0;
const Vector4d obj_reflection_color(0.7, 0.7, 0.7, 0);
const Vector4d obj_refraction_color(0.7, 0.7, 0.7, 0);

// Refraction parameters
bool enable_refraction = false;
double src_refr_index = 1; // source refractive index (vacuum/air assumed)
double obj_refr_index = 1.5; // sphere refractive index (glass value used)
double refr_critical_angle = 1.57; // vacuum to glass

// Precomputed (or otherwise) gradient vectors at each grid node
const int grid_size = 20;
std::vector<std::vector<Vector2d>> grid;

//Lights
std::vector<Vector3d> light_positions;
std::vector<Vector4d> light_colors;
//Ambient light
const Vector4d ambient_light(0.2, 0.2, 0.2, 0);

// Depth-of-field parameters
// Set dof_sample_radius = 0.0 to turn off dof
int dof_n_samples = 10;
double dof_sample_radius = 0.06;
std::vector<Vector3d> dof_noise;

//Fills the different arrays
void setup_scene()
{
    grid.resize(grid_size + 1);
    for (int i = 0; i < grid_size + 1; ++i)
    {
        grid[i].resize(grid_size + 1);
        for (int j = 0; j < grid_size + 1; ++j)
            grid[i][j] = Vector2d::Random().normalized();
    }

    // DOF random points in dof_sample_radius around the camera
    std::uniform_real_distribution<double> unif(0.0,dof_sample_radius);
    std::default_random_engine re;
    for (int i=0; i<dof_n_samples; i++) {
        double rand_radius = unif(re);
        Vector2d rand_point = Vector2d::Random().normalized() * rand_radius;
        double noisy_x = camera_position.x() + rand_point.x();
        double noisy_y = camera_position.y() + rand_point.y();
        dof_noise.emplace_back(Vector3d(noisy_x, noisy_y, 0));
    }


    //Spheres
    sphere_centers.emplace_back(10, 0, 1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(7, 0.05, -1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(4, 0.1, 1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(1, 0.2, -1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(-2, 0.4, 1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(-5, 0.8, -1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(-8, 1.6, 1);
    sphere_radii.emplace_back(1);

    //parallelograms
    parallelograms.emplace_back();
    parallelograms.back() << -100, 100, -100,
        -1.25, 0, -1.2,
        -100, -100, 100;

    //Lights
    light_positions.emplace_back(8, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(6, -8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(4, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(2, -8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(0, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(-2, -8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(-4, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);
}

//We need to make this function visible
Vector4d shoot_ray(const Vector3d &ray_origin, const Vector3d &ray_direction, int max_bounce);

////////////////////////////////////////////////////////////////////////////////
// Perlin noise code
////////////////////////////////////////////////////////////////////////////////

// Function to linearly interpolate between a0 and a1
// Weight w should be in the range [0.0, 1.0]
double lerp(double a0, double a1, double w)
{
    assert(w >= 0-0.001); // Small tolerance for floating point errors
    assert(w <= 1+0.001);

    return a0 + w * (a1-a0); // Linear interpolation
//    return (a1 - a0) * (3.0 - w * 2.0) * w * w + a0; // Cubic Interpolation
}

// Computes the dot product of the distance and gradient vectors.
double dotGridGradient(int ix, int iy, double x, double y)
{
//    extern std::vector<std::vector<Vector2d>> grid;
    double dx = x - (double)ix;
    double dy = y - (double)iy;
    return dx*grid[ix][iy][0] + dy*grid[ix][iy][1];
}

// Compute Perlin noise at coordinates x, y
double perlin(double x, double y)
{
    int x0 = int(x);
    int x1 = x0 + 1;
    int y0 = int(y);
    int y1 = y0 + 1;

    // Determine interpolation weights
    double sx = x - x0;
    double sy = y - y0;

    // Interpolate between grid point gradients
    double n0 = dotGridGradient(x0, y0, x, y);
    double n1 = dotGridGradient(x1, y0, x, y);

    double ix0 = lerp(n0, n1, sx);

    n0 = dotGridGradient(x0, y1, x, y);
    n1 = dotGridGradient(x1, y1, x, y);

    double ix1 = lerp(n0, n1, sx);
    double value = lerp(ix0, ix1, sy);

    return value;
}

Vector4d procedural_texture(const double tu, const double tv)
{
    assert(tu >= 0-0.001); // Small tolerance for floating point errors
    assert(tv >= 0-0.001);

    assert(tu <= 1+0.001); // Small tolerance for floating point errors
    assert(tv <= 1+0.001);

    const double color = (perlin(tu * grid_size, tv * grid_size) + 1) / 2;
    return Vector4d(0, color, 0, 0);

    //Example for checkerboard texture
//    const double color = (int(tu * grid_size) + int(tv * grid_size)) % 2 == 0 ? 0 : 1;
//    return Vector4d(0, color, 0, 0);
}

////////////////////////////////////////////////////////////////////////////////
// Intersection code
////////////////////////////////////////////////////////////////////////////////

//Compute the intersection between a ray and a sphere, return -1 if no intersection
double ray_sphere_intersection(const Vector3d &ray_origin, const Vector3d &ray_direction, int index, Vector3d &p, Vector3d &N)
{
    //return t or -1 if no intersection

    const Vector3d sphere_center = sphere_centers[index];
    const double sphere_radius = sphere_radii[index];

    double t = -1;

    double A = ray_direction.dot(ray_direction);
    double B = 2 * ray_direction.dot(ray_origin-sphere_center);
    double C = (ray_origin-sphere_center).dot(ray_origin-sphere_center) - sphere_radius*sphere_radius;
    double discriminant = B*B - 4*A*C;
    if (discriminant < 0)
    {
        return -1;
    }
    else
    {
        t = std::min((-B-std::sqrt(discriminant))/(2*A), (-B+std::sqrt(discriminant))/(2*A));
        p = ray_origin + t*ray_direction;
        N = (p-sphere_center).normalized();
        return t;
    }

    return -1;
}

//Compute the intersection between a ray and a paralleogram, return -1 if no intersection
double ray_parallelogram_intersection(const Vector3d &ray_origin, const Vector3d &ray_direction, int index, Vector3d &p, Vector3d &N)
{
    //return t or -1 if no intersection

    const Vector3d pgram_origin = parallelograms[index].col(0);
    const Vector3d pgram_u = parallelograms[index].col(1) - pgram_origin;
    const Vector3d pgram_v = parallelograms[index].col(2) - pgram_origin;

    Matrix3d A;
    Vector3d b;
    A << pgram_u.x(), pgram_v.x(), -ray_direction.x(), pgram_u.y(), pgram_v.y(), -ray_direction.y(), pgram_u.z(), pgram_v.z(), -ray_direction.z();
    b = ray_origin - pgram_origin;

    Vector3d x = A.colPivHouseholderQr().solve(b);
    double p1 = x.x();
    double q1 = x.y();
    double t = x.z();
    if (!(p1>=0 && p1<=1 && q1>=0 && q1<=1 && t>0))
    {
        return -1;
    }

    p = ray_origin + t*ray_direction;
    N = pgram_u.cross(pgram_v).normalized();
    if (N.dot(-1*ray_direction) < 0) {
        N = -1 * N;
    }
    return t;
}

//Finds the closest intersecting object returns its index
//In case of intersection it writes into p and N (intersection point and normals)
int find_nearest_object(const Vector3d &ray_origin, const Vector3d &ray_direction, Vector3d &p, Vector3d &N)
{
    // Find the object in the scene that intersects the ray first
    // we store the index and the 'closest_t' to their expected values
    int closest_index = -1;
    double closest_t = std::numeric_limits<double>::max(); //closest t is "+ infinity"

    Vector3d tmp_p, tmp_N;
    for (int i = 0; i < sphere_centers.size(); ++i)
    {
        //returns t and writes on tmp_p and tmp_N
        const double t = ray_sphere_intersection(ray_origin, ray_direction, i, tmp_p, tmp_N);
        //We have intersection
        if (t >= 0)
        {
            //The point is before our current closest t
            if (t < closest_t)
            {
                closest_index = i;
                closest_t = t;
                p = tmp_p;
                N = tmp_N;
            }
        }
    }

    for (int i = 0; i < parallelograms.size(); ++i)
    {
        //returns t and writes on tmp_p and tmp_N
        const double t = ray_parallelogram_intersection(ray_origin, ray_direction, i, tmp_p, tmp_N);
        //We have intersection
        if (t >= 0)
        {
            //The point is before our current closest t
            if (t < closest_t)
            {
                closest_index = sphere_centers.size() + i;
                closest_t = t;
                p = tmp_p;
                N = tmp_N;
            }
        }
    }

    return closest_index;
}

////////////////////////////////////////////////////////////////////////////////
// Raytracer code
////////////////////////////////////////////////////////////////////////////////

//Checks if the light is visible
bool is_light_visible(const Vector3d &ray_origin, const Vector3d &ray_direction, const Vector3d &light_position)
{
    Vector3d p, N;
    const int nearest_object = find_nearest_object(ray_origin + epsilon*ray_direction, ray_direction, p, N);
    if (nearest_object < 0) {
        return true;
    }
    return (light_position-ray_origin).norm() < (p-ray_origin).norm();
}

Vector4d shoot_ray(const Vector3d &ray_origin, const Vector3d &ray_direction, int max_bounce)
{
    //Intersection point and normal, these are output of find_nearest_object
    Vector3d p, N;

    const int nearest_object = find_nearest_object(ray_origin, ray_direction, p, N);

    if (nearest_object < 0)
    {
        // Return a transparent color
        return Vector4d(0, 0, 0, 0);
    }

    // Ambient light contribution
    const Vector4d ambient_color = obj_ambient_color.array() * ambient_light.array();

    // Punctual lights contribution (direct lighting)
    Vector4d lights_color(0, 0, 0, 0);
    for (int i = 0; i < light_positions.size(); ++i)
    {
        const Vector3d &light_position = light_positions[i];
        const Vector4d &light_color = light_colors[i];

        const Vector3d Li = (light_position - p).normalized();

        if (!is_light_visible(p, Li, light_position)) {
            continue;
        }

        Vector4d diff_color = obj_diffuse_color;

        if (nearest_object == 4)
        {
            //Compute UV coodinates for the point on the sphere
            const double x = p(0) - sphere_centers[nearest_object][0];
            const double y = p(1) - sphere_centers[nearest_object][1];
            const double z = p(2) - sphere_centers[nearest_object][2];
            const double tu = acos(z / sphere_radii[nearest_object]) / 3.1415;
            const double tv = (3.1415 + atan2(y, x)) / (2 * 3.1415);

            diff_color = procedural_texture(tu, tv);
        }

        // Diffuse contribution
        const Vector4d diffuse = diff_color * std::max(Li.dot(N), 0.0);

        const Vector3d h = ( Li + (camera_position-p).normalized() ).normalized();
        // Specular contribution, use obj_specular_color
        const Vector4d specular = obj_specular_color * pow(std::max(0.0, N.dot(h)), obj_specular_exponent);

        // Attenuate lights according to the squared distance to the lights
        const Vector3d D = light_position - p;
        lights_color += (diffuse + specular).cwiseProduct(light_color) / D.squaredNorm();
    }

    Vector4d refl_color = obj_reflection_color;
    if (nearest_object == 4)
    {
        refl_color = Vector4d(0.5, 0.5, 0.5, 0);
    }

    Vector4d reflection_color(0, 0, 0, 0);
    if (max_bounce>0) {
        Vector3d new_ray_direction = ray_direction - 2*(ray_direction.dot(N))*N;
        reflection_color = shoot_ray(p + epsilon*new_ray_direction, new_ray_direction, max_bounce-1).cwiseProduct(refl_color);
    }

    Vector4d refraction_color(0, 0, 0, 0);
    if (enable_refraction) {
        // Compute the color of the refracted ray and add its contribution to the current point color.
        double cos_in = N.dot((-ray_direction).normalized()); // cos of angle of incidence
        if (acos(cos_in) < refr_critical_angle) {
            // cos of angle of refraction
            double cos_refr = std::sqrt(1.0-(std::pow(src_refr_index,2)*(1-std::pow(cos_in,2))/std::pow(obj_refr_index,2)));
            Vector3d new_ray_direction = ((src_refr_index/obj_refr_index) * (ray_direction + N*cos_in)) - N*cos_refr;
            refraction_color = shoot_ray(p + epsilon * new_ray_direction, new_ray_direction,
                                                 max_bounce - 1).cwiseProduct(obj_refraction_color);
        }
    }

    // Rendering equation
    Vector4d C = ambient_color + lights_color + reflection_color + refraction_color;

    //Set alpha to 1
    C(3) = 1;

    return C;
}

////////////////////////////////////////////////////////////////////////////////

void raytrace_scene()
{
    std::cout << "Simple ray tracer." << std::endl;

    int w = 800;
    int h = 400;
    MatrixXd R = MatrixXd::Zero(w, h);
    MatrixXd G = MatrixXd::Zero(w, h);
    MatrixXd B = MatrixXd::Zero(w, h);
    MatrixXd A = MatrixXd::Zero(w, h); // Store the alpha mask

    // The camera always points in the direction -z
    // The sensor grid is at a distance 'focal_length' from the camera center,
    // and covers an viewing angle given by 'field_of_view'.
    double aspect_ratio = double(w) / double(h);
    double image_y = focal_length * tan(field_of_view/2);
    double image_x = aspect_ratio * image_y;

    // The pixel grid through which we shoot rays is at a distance 'focal_length'
    const Vector3d image_origin(-image_x, image_y, -image_z);
    const Vector3d x_displacement(2.0 / w * image_x, 0, 0);
    const Vector3d y_displacement(0, -2.0 / h * image_y, 0);

    for (unsigned i = 0; i < w; ++i)
    {
        for (unsigned j = 0; j < h; ++j)
        {
            const Vector3d pixel_center = image_origin + (i + 0.5) * x_displacement + (j + 0.5) * y_displacement;

            // Prepare the ray
            Vector3d ray_origin;
            Vector3d ray_direction;

            // Depth-of-field
            for (int k=0; k<dof_n_samples; k++) {
                if (is_perspective)
                {
                    ray_origin = camera_position + dof_noise[k];
                    ray_direction = (pixel_center-ray_origin).normalized();
                }
                else
                {
                    // Orthographic camera
                    ray_origin = camera_position + Vector3d(pixel_center[0], pixel_center[1], 0) + dof_noise[k];
                    ray_direction = Vector3d(0, 0, -1);
                }

                const Vector4d C = shoot_ray(ray_origin, ray_direction, max_bounce);
                R(i, j) += C(0);
                G(i, j) += C(1);
                B(i, j) += C(2);
                A(i, j) += C(3);
            }
            R(i, j) /= dof_n_samples;
            G(i, j) /= dof_n_samples;
            B(i, j) /= dof_n_samples;
            A(i, j) /= dof_n_samples;
        }
    }

    // Save to png
    write_matrix_to_png(R, G, B, A, filename);
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    setup_scene();

    raytrace_scene();
    return 0;
}
