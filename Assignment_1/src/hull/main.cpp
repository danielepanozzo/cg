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

struct Compare {
	Point p0; // Leftmost point of the poly
	bool operator ()(const Point &p1, const Point &p2) {
        double cross_product = det(p1-p0, p2-p0);
        if(cross_product==0) {
            return norm(p2-p0) >= norm(p1-p0);
        }
		return cross_product>0;
	}
};

bool inline salientAngle(Point &a, Point &b, Point &c) {
	return det(b-a, c-b)>=0;
}

////////////////////////////////////////////////////////////////////////////////

Polygon convex_hull(std::vector<Point> &points) {
	Compare order;

    // Scanning all points to find the point p0 with lowest y coordinate
    double p0_y = points.at(0).imag();
    double p0_x = points.at(0).real();
    for (int i=0; i<points.size(); i++) {
        if(points.at(i).imag() < p0_y) {
            p0_y = points.at(i).imag();
            p0_x = points.at(i).real();
        } else if (points.at(i).imag() == p0_y && points.at(i).real() < p0_x) {
            p0_y = points.at(i).imag();
            p0_x = points.at(i).real();
        }
    }
    order.p0 = Point(p0_x, p0_y);

	std::sort(points.begin(), points.end(), order);
	Polygon hull;
	for(int i=0; i<points.size(); i++) {
        while(hull.size()>2 && !salientAngle(hull.at(hull.size()-2), hull.back(), points.at(i))) {
            hull.pop_back();
        }
        hull.push_back(points.at(i));
    }

	return hull;
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

void save_obj(const std::string &filename, Polygon &poly) {
	std::ofstream out(filename);
	if (!out.is_open()) {
		throw std::runtime_error("failed to open file " + filename);
	}
	out << std::fixed;
	for (const auto &v : poly) {
		out << "v " << v.real() << ' ' << v.imag() << " 0\n";
	}
	for (size_t i = 0; i < poly.size(); ++i) {
		out << "l " << i+1 << ' ' << 1+(i+1)%poly.size() << "\n";
	}
	out << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char * argv[]) {
	if (argc <= 2) {
		std::cerr << "Usage: " << argv[0] << " points.xyz output.obj" << std::endl;
	}
	std::vector<Point> points = load_xyz(argv[1]);
	Polygon hull = convex_hull(points);
	save_obj(argv[2], hull);
	return 0;
}
