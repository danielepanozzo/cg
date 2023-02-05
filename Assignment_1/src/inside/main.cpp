////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

typedef std::complex<double> Point;
typedef std::vector<Point> Polygon;

double inline det(const Point &u, const Point &v) {
    return u.real()*v.imag() - u.imag()*v.real();
}

// Return true iff [a,b] intersects [c,d], and store the intersection in ans
bool intersect_segment(const Point &a, const Point &b, const Point &c, const Point &d, Point &ans) {
    // representing a line as y = mx + c
	double m1 = (b.imag()-a.imag()) / (b.real()-a.real());
    double m2 = (d.imag()-c.imag()) / (d.real()-c.real());
    double c1 = (b.real()*a.imag() - a.real()*b.imag()) / (b.real()-a.real());
    double c2 = (d.real()*c.imag() - c.real()*d.imag()) / (d.real()-c.real());
    if(m1==m2) {
        return false;
    }
    double x_int = (c2-c1)/(m1-m2);
    double y_int = (m1*c2-m2*c1)/(m1-m2);
    if( (x_int > std::min(a.real(), b.real())) &&
        (x_int < std::max(a.real(), b.real())) &&
        (x_int > std::min(c.real(), d.real())) &&
        (x_int < std::max(c.real(), d.real())) ) {
        ans = Point(x_int, y_int);
        return true;
    }
	return false;
}

////////////////////////////////////////////////////////////////////////////////

bool is_inside(const Polygon &poly, const Point &query) {
	// 1. Compute bounding box and set coordinate of a point outside the polygon
	double max_x = poly.at(0).real();
    double max_y = poly.at(0).imag();
    for(auto v : poly) {
        max_x = std::max(max_x, v.real());
        max_y = std::max(max_y, v.imag());
    }
	Point outside(2*max_x, 2*max_y);

	// 2. Cast a ray from the query point to the 'outside' point, count number of intersections
	int n_intersections = 0;
    Point ans;
    for (int i=0; i<poly.size(); i++) {
        if(intersect_segment(query, outside, poly.at(i), poly.at((i+1)%poly.size()), ans)) {
            n_intersections += 1;
        }
    }
	return n_intersections % 2 == 1;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Point> load_xyz(const std::string &filename) {
    std::vector<Point> points;
    std::ifstream in(filename);
    std::string line;
    std::getline(in, line);
    int n_points = std::stoi(line);
    float x, y, z;
    for (int i=0; i<n_points; i++) {
        in >> x >> y >> z;
        points.push_back(Point(x, y));
    }
    return points;
}

Polygon load_obj(const std::string &filename) {
    std::vector<Point> input_polygon;
	std::ifstream in(filename);
    std::string line;
    float x, y, z;
    std::string v;
    while(getline(in, line)) {
        if(line[0]=='v') {
            std::istringstream iss( line );
            iss >> v >> x >> y >> z;
        }
        input_polygon.push_back(Point(x, y));
    }
	return input_polygon;
}

void save_xyz(const std::string &filename, const std::vector<Point> &points) {
	std::ofstream out(filename);
    if (!out.is_open()) {
        throw std::runtime_error("failed to open file " + filename);
    }
    out << std::fixed;
    out << points.size() << "\n";
    for (const auto &v : points) {
        out << v.real() << ' ' << v.imag() << " 0\n";
    }
    out << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char * argv[]) {
	if (argc <= 3) {
		std::cerr << "Usage: " << argv[0] << " points.xyz poly.obj result.xyz" << std::endl;
	}
	std::vector<Point> points = load_xyz(argv[1]);
	Polygon poly = load_obj(argv[2]);
	std::vector<Point> result;
	for (size_t i = 0; i < points.size(); ++i) {
		if (is_inside(poly, points[i])) {
			result.push_back(points[i]);
		}
	}
	save_xyz(argv[3], result);
	return 0;
}
