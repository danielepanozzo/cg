////////////////////////////////////////////////////////////////////////////////
// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <float.h>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

////////////////////////////////////////////////////////////////////////////////
// Class to store tree
////////////////////////////////////////////////////////////////////////////////
class AABBTree
{
public:
    class Node
    {
    public:
        AlignedBox3d bbox;
        int parent;   // Index of the parent node (-1 for root)
        int left;     // Index of the left child (-1 for a leaf)
        int right;    // Index of the right child (-1 for a leaf)
        int triangle; // Index of the node triangle (-1 for internal nodes)
    };

    std::vector<Node> nodes;
    int root;

    AABBTree() = default;                           // Default empty constructor
    AABBTree(const MatrixXd &V, const MatrixXi &F); // Build a BVH from an existing mesh
    int build_bvh_top_down(const MatrixXd &V, const MatrixXi &F, const MatrixXd &centroids, int* s, int* e);
};

////////////////////////////////////////////////////////////////////////////////
// Scene setup, global variables
////////////////////////////////////////////////////////////////////////////////
const std::string data_dir = DATA_DIR;
const std::string filename("raytrace.png");
const std::string mesh_filename(data_dir + "dodeca.off");

//Camera settings
const double focal_length = 2;
const double field_of_view = 0.7854; //45 degrees
const bool is_perspective = true;
const Vector3d camera_position(0, 0, 2);

// Triangle Mesh
MatrixXd vertices; // n x 3 matrix (n points)
MatrixXi facets;   // m x 3 matrix (m triangles)
AABBTree bvh;

//Material for the object, same material for all objects
const Vector4d obj_ambient_color(0.0, 0.5, 0.0, 0);
const Vector4d obj_diffuse_color(0.5, 0.5, 0.5, 0);
const Vector4d obj_specular_color(0.2, 0.2, 0.2, 0);
const double obj_specular_exponent = 256.0;
const Vector4d obj_reflection_color(0.7, 0.7, 0.7, 0);

// Precomputed (or otherwise) gradient vectors at each grid node
const int grid_size = 20;
std::vector<std::vector<Vector2d>> grid;

//Lights
std::vector<Vector3d> light_positions;
std::vector<Vector4d> light_colors;
//Ambient light
const Vector4d ambient_light(0.2, 0.2, 0.2, 0);

//Fills the different arrays
void setup_scene()
{
    //Loads file
    std::ifstream in(mesh_filename);
    std::string token;
    in >> token;
    int nv, nf, ne;
    in >> nv >> nf >> ne;
    vertices.resize(nv, 3);
    facets.resize(nf, 3);
    for (int i = 0; i < nv; ++i)
    {
        in >> vertices(i, 0) >> vertices(i, 1) >> vertices(i, 2);
    }
    for (int i = 0; i < nf; ++i)
    {
        int s;
        in >> s >> facets(i, 0) >> facets(i, 1) >> facets(i, 2);
        assert(s == 3);
    }

    //setup tree
    bvh = AABBTree(vertices, facets);

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

////////////////////////////////////////////////////////////////////////////////
// BVH Code
////////////////////////////////////////////////////////////////////////////////

AlignedBox3d bbox_from_triangle(const Vector3d &a, const Vector3d &b, const Vector3d &c)
{
    AlignedBox3d box;
    box.extend(a);
    box.extend(b);
    box.extend(c);
    return box;
}

AABBTree::AABBTree(const MatrixXd &V, const MatrixXi &F)
{
    // Compute the centroids of all the triangles in the input mesh
    MatrixXd centroids(F.rows(), V.cols());
    centroids.setZero();
    for (int i = 0; i < F.rows(); ++i)
    {
        for (int k = 0; k < F.cols(); ++k)
        {
            centroids.row(i) += V.row(F(i, k));
        }
        centroids.row(i) /= F.cols();
    }

    // Split each set of primitives into 2 sets of roughly equal size,
    // based on sorting the centroids along one direction or another.
    int triangle_idx[F.rows()];
    for (int i=0; i<F.rows(); i++) {
        triangle_idx[i] = i;
    }
    root = build_bvh_top_down(V, F, centroids, triangle_idx, triangle_idx+F.rows());
}

int AABBTree::build_bvh_top_down(const MatrixXd &V, const MatrixXi &F, const MatrixXd &centroids, int* s, int* e) {
    Node new_node;
    new_node.parent = -1;


    if ((e - s) == 1) {
        // This indicates a leaf node
        new_node.left = -1;
        new_node.right = -1;
        new_node.triangle = *s;

        new_node.bbox = bbox_from_triangle(V.row(F(*s, 0)), V.row(F(*s, 1)), V.row(F(*s, 2)));
        nodes.emplace_back(new_node);
        return nodes.size()-1;
    }
    // Decide axis to sort the centroids
    AlignedBox3d b;
    for (int *i = s; i < e; i++) {
        b.extend(centroids.row(*i).transpose());
    }
    Vector3d b_sizes = b.sizes();
    int split_axis = 0;
    double curr_max = b_sizes[0];
    for (int i = 1; i < centroids.cols(); i++) {
        if (b_sizes[i] > curr_max) {
            split_axis = i;
            curr_max = b_sizes[i];
        }
    }

    // Sort based on split_axis and divide points into two sets
    // TODO: enable sorting
//    std::sort(s, e, [&](int a, int b) { return centroids(a, split_axis) < centroids(b, split_axis); });
    int *m = s + (e - s + 1) / 2;
    new_node.left = build_bvh_top_down(V, F, centroids, s, m);
    new_node.right = build_bvh_top_down(V, F, centroids, m, e);
    new_node.bbox = AlignedBox3d();
    new_node.bbox.extend(nodes[new_node.left].bbox);
    new_node.bbox.extend(nodes[new_node.right].bbox);
    new_node.triangle = -1;
    nodes.emplace_back(new_node);
    nodes[new_node.left].parent = nodes.size()-1;
    nodes[new_node.right].parent = nodes.size()-1;

    return nodes.size()-1;
}

////////////////////////////////////////////////////////////////////////////////
// Intersection code
////////////////////////////////////////////////////////////////////////////////

double ray_triangle_intersection(const Vector3d &ray_origin, const Vector3d &ray_direction, const Vector3d &a, const Vector3d &b, const Vector3d &c, Vector3d &p, Vector3d &N)
{
    // Compute whether the ray intersects the given triangle.
    // If you have done the parallelogram case, this should be very similar to it.


    Matrix3d A_eq;
    Vector3d b_eq;
    A_eq << a.x()-b.x(), a.x()-c.x(), ray_direction.x(),
            a.y()-b.y(), a.y()-c.y(), ray_direction.y(),
            a.z()-b.z(), a.z()-c.z(), ray_direction.z();
    b_eq = a - ray_origin;

    Vector3d x = A_eq.colPivHouseholderQr().solve(b_eq);
    double u = x.x();
    double v = x.y();
    double t = x.z();
    if (!(u>=0 && v>=0 && (u+v)<=1 && t>0)) {
        return -1;
    }

    p = ray_origin + t*ray_direction;
    N = (b-a).cross(c-a).normalized();
    if (N.dot(-1*ray_direction) < 0) {
        N = -1 * N;
    }
    return t;
}

bool ray_box_intersection(const Vector3d &ray_origin, const Vector3d &ray_direction, const AlignedBox3d &box)
{
    // Compute whether the ray intersects the given box.
    // we are not testing with the real surface here anyway.

    // Reference: Fundamentals of Computer Graphics, 4th edition, Section 12.3.1
    double t_xmin, t_xmax, t_ymin, t_ymax, t_zmin, t_zmax;
    Vector3d inv_ray_dir;
    inv_ray_dir.x() = 1.0/ray_direction.x();
    inv_ray_dir.y() = 1.0/ray_direction.y();
    inv_ray_dir.z() = 1.0/ray_direction.z();

    double temp;
    t_xmin = (box.min().x() - ray_origin.x()) * inv_ray_dir.x();
    t_xmax = (box.max().x() - ray_origin.x()) * inv_ray_dir.x();
    if (t_xmax < t_xmin) {
        temp = t_xmin;
        t_xmin = t_xmax;
        t_xmax = temp;
    }

    t_ymin = (box.min().y() - ray_origin.y()) * inv_ray_dir.y();
    t_ymax = (box.max().y() - ray_origin.y()) * inv_ray_dir.y();
    if (t_ymax < t_ymin) {
        temp = t_ymin;
        t_ymin = t_ymax;
        t_ymax = temp;
    }

    t_zmin = (box.min().z() - ray_origin.z()) * inv_ray_dir.z();
    t_zmax = (box.max().z() - ray_origin.z()) * inv_ray_dir.z();
    if (t_zmax < t_zmin) {
        temp = t_zmin;
        t_zmin = t_zmax;
        t_zmax = temp;
    }


//    if(inv_ray_dir.x() >= 0) {
//        t_xmin = (box.min().x() - ray_origin.x()) * inv_ray_dir.x();
//        t_xmax = (box.max().x() - ray_origin.x()) * inv_ray_dir.x();
//    } else {
//        t_xmin = (box.max().x() - ray_origin.x()) * inv_ray_dir.x();
//        t_xmax = (box.min().x() - ray_origin.x()) * inv_ray_dir.x();
//    }
//    if(inv_ray_dir.y() >= 0) {
//        t_ymin = (box.min().y() - ray_origin.y()) * inv_ray_dir.y();
//        t_ymax = (box.max().y() - ray_origin.y()) * inv_ray_dir.y();
//    } else {
//        t_ymin = (box.max().y() - ray_origin.y()) * inv_ray_dir.y();
//        t_ymax = (box.min().y() - ray_origin.y()) * inv_ray_dir.y();
//    }
//    if(inv_ray_dir.z() >= 0) {
//        t_zmin = (box.min().z() - ray_origin.z()) * inv_ray_dir.z();
//        t_zmax = (box.max().z() - ray_origin.z()) * inv_ray_dir.z();
//    } else {
//        t_zmin = (box.max().z() - ray_origin.z()) * inv_ray_dir.z();
//        t_zmax = (box.min().z() - ray_origin.z()) * inv_ray_dir.z();
//    }
    double t_min = std::max(t_xmin, std::max(t_ymin, t_zmin));
    double t_max = std::min(t_xmax, std::min(t_ymax, t_zmax));
    if (t_max < 0) {
        return false;
    }
    return (t_max >= t_min);
}

bool find_nearest_object_brute_force(const Vector3d &ray_origin, const Vector3d &ray_direction, Vector3d &p, Vector3d &N) {
    double t_min = DBL_MAX;
    Vector3d nearest_triangle_p;
    Vector3d nearest_traingle_N;

    double t;
    Vector3d p1;
    Vector3d N1;
    for (int i=0; i<facets.rows(); i++) {
        t = ray_triangle_intersection(
                ray_origin, ray_direction,
                vertices.row(facets(i, 0)),
                vertices.row(facets(i, 1)),
                vertices.row(facets(i, 2)),
                p1, N1
        );
        if(t>=0 && t<t_min) {
            t_min = t;
            nearest_triangle_p = p1;
            nearest_traingle_N = N1;
        }
    }
    if (t_min == DBL_MAX) {
        return false;
    }
    p = nearest_triangle_p;
    N = nearest_traingle_N;
    return true;
}

bool find_nearest_object_bvh(const int curr_node_id, const Vector3d &ray_origin, const Vector3d &ray_direction, Vector3d &p, Vector3d &N) {
    if (bvh.nodes[curr_node_id].left == -1 && bvh.nodes[curr_node_id].right == -1) {
        // Leaf node, check intersection with triangle
        return ray_triangle_intersection(
                ray_origin, ray_direction,
                vertices.row(facets(bvh.nodes[curr_node_id].triangle, 0)),
                vertices.row(facets(bvh.nodes[curr_node_id].triangle, 1)),
                vertices.row(facets(bvh.nodes[curr_node_id].triangle, 2)),
                p, N);
    }
    if (ray_box_intersection(ray_origin, ray_direction, bvh.nodes[curr_node_id].bbox)) {
        Vector3d* p1 = new Vector3d();
        Vector3d* p2 = new Vector3d();
        Vector3d* N1 = new Vector3d();
        Vector3d* N2 = new Vector3d();

        bool found_left = find_nearest_object_bvh(bvh.nodes[curr_node_id].left, ray_origin, ray_direction, *p1, *N1);
        bool found_right = find_nearest_object_bvh(bvh.nodes[curr_node_id].right, ray_origin, ray_direction, *p2, *N2);

        if (found_left && found_right) {
            // pick nearest
            if ((*p1-ray_origin).norm() >= (*p2-ray_origin).norm()) {
                p = *p2;
                N = *N2;
                delete p1;
                delete N1;
            } else {
                p = *p1;
                N = *N1;
                delete p2;
                delete N2;
            }
            return true;
        } else if (found_left) {
            p = *p1;
            N = *N1;
            delete p2;
            delete N2;
            return true;
        } else if (found_right) {
            p = *p2;
            N = *N2;
            delete p1;
            delete N1;
            return true;
        } else {
            delete p1;
            delete N1;
            delete p2;
            delete N2;
            return false;
        }
    } else {
        return false;
    }
}

//Finds the closest intersecting object returns its index
//In case of intersection it writes into p and N (intersection point and normals)
bool find_nearest_object(const Vector3d &ray_origin, const Vector3d &ray_direction, Vector3d &p, Vector3d &N)
{
    bool use_bvh = true;

    if (!use_bvh) {
        // Method (1) Traverse every triangle and return the closest hit.
        return find_nearest_object_brute_force(ray_origin, ray_direction, p, N);
    } else {
        // Method (2): Traverse the BVH tree and test the intersection with a
        // triangles at the leaf nodes that intersects the input ray.
        return find_nearest_object_bvh(bvh.root, ray_origin, ray_direction, p, N);
    }
}

////////////////////////////////////////////////////////////////////////////////
// Raytracer code
////////////////////////////////////////////////////////////////////////////////

Vector4d shoot_ray(const Vector3d &ray_origin, const Vector3d &ray_direction)
{
    //Intersection point and normal, these are output of find_nearest_object
    Vector3d p, N;

    const bool nearest_object = find_nearest_object(ray_origin, ray_direction, p, N);

    if (!nearest_object)
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

        Vector4d diff_color = obj_diffuse_color;

        // Diffuse contribution
        const Vector3d Li = (light_position - p).normalized();
        const Vector4d diffuse = diff_color * std::max(Li.dot(N), 0.0);

        // Specular contribution
        const Vector3d Hi = (Li - ray_direction).normalized();
        const Vector4d specular = obj_specular_color * std::pow(std::max(N.dot(Hi), 0.0), obj_specular_exponent);
        // Vector3d specular(0, 0, 0);

        // Attenuate lights according to the squared distance to the lights
        const Vector3d D = light_position - p;
        lights_color += (diffuse + specular).cwiseProduct(light_color) / D.squaredNorm();
    }

    // Rendering equation
    Vector4d C = ambient_color + lights_color;

    //Set alpha to 1
    C(3) = 1;

    return C;
}

////////////////////////////////////////////////////////////////////////////////

void raytrace_scene()
{
    std::cout << "Simple ray tracer." << std::endl;

    int w = 640;
    int h = 480;
    MatrixXd R = MatrixXd::Zero(w, h);
    MatrixXd G = MatrixXd::Zero(w, h);
    MatrixXd B = MatrixXd::Zero(w, h);
    MatrixXd A = MatrixXd::Zero(w, h); // Store the alpha mask

    // The camera always points in the direction -z
    // The sensor grid is at a distance 'focal_length' from the camera center,
    // and covers an viewing angle given by 'field_of_view'.
    double aspect_ratio = double(w) / double(h);
    double image_y = focal_length * std::tan(field_of_view/2);
    double image_x = aspect_ratio * image_y;

    // The pixel grid through which we shoot rays is at a distance 'focal_length'
    const Vector3d image_origin(-image_x, image_y, camera_position[2] - focal_length);
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

            if (is_perspective)
            {
                // Perspective camera
                ray_origin = camera_position;
                ray_direction = (pixel_center - camera_position).normalized();
            }
            else
            {
                // Orthographic camera
                ray_origin = pixel_center;
                ray_direction = Vector3d(0, 0, -1);
            }

            const Vector4d C = shoot_ray(ray_origin, ray_direction);
            R(i, j) = C(0);
            G(i, j) = C(1);
            B(i, j) = C(2);
            A(i, j) = C(3);
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