#include <vector>
#include <cmath>
#include <math.h>
#include <random>
#include <algorithm>

#include "objects.h"

using namespace std;

// Vector definition
Vector::Vector(double x, double y, double z) {
    coords[0] = x;
    coords[1] = y;
    coords[2] = z;
};

Vector& Vector::operator+=(const Vector& b) {
    coords[0] += b[0];
    coords[1] += b[1];
    coords[2] += b[2];
    return *this;
}

const double& Vector::operator[](int i) const {return coords[i]; }
double& Vector::operator[](int i) { return coords[i]; }


Vector operator+(const Vector& a, const Vector &b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector operator-(const Vector& a, const Vector &b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector operator*(double a, const Vector& v) {
    return  Vector(a * v[0], a * v[1], a * v[2]);
}

Vector operator*(const Vector& v, double a) {
    return Vector(a * v[0], a * v[1], a * v[2]);
}

Vector operator*(const Vector& v, const Vector& a) {      // tensor product
    return Vector(a[0] * v[0], a[1] * v[1], a[2] * v[2]);
}

Vector operator/(const Vector& v, double a) {
    return Vector(v[0] / a, v[1] / a, v[2] / a); 
}

Vector scale(double a, const Vector& v) {
    return  Vector(a * v[0], a * v[1], a * v[2]);
}

double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double norm(const Vector& a) {
    return sqrt(dot(a, a));
}

Vector unit(const Vector& a) {
    return scale(1 / norm(a), a);
}

Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

// Polygon 

Vector Polygon::center() {

    Vector center(0, 0, 0);

    for (auto &vert : vertices)
        center += vert;

    return center / vertices.size(); 

}

bool is_inside(Vector const &P, Vector const &u, Vector const &v, Vector const &inner_point) {

    Vector N = Vector(v[1] - u[1], u[0] - v[0], 0);

    // makes sure N points outwards
    // N = (dot(u - inner_point, N) >= 0) ? N : -1 * N;

    return (dot(P - u, N) <= 0) ? true : false;
}

Vector intersect(Vector const &A, Vector const &B, Vector const &u, Vector const &v) {

    Vector N = Vector(v[1] - u[1], u[0] - v[0], 0);
    double t = dot(u - A, N) / dot(B - A, N);

    return A + t * (B - A);

}

Polygon clip_polygon(Polygon subject_polygon, Polygon &clip_polygon) {

    int clip_n = clip_polygon.vertices.size();

    Vector clip_center = clip_polygon.center();

    Polygon out_polygon;

    for (int clip_i = 0; clip_i < clip_n; ++clip_i) {

        out_polygon = Polygon();
        int subj_n = subject_polygon.vertices.size();

        for (int subj_i = 0; subj_i < subj_n; ++subj_i) {

            Vector clip_vert1 = clip_polygon.vertices[clip_i];
            Vector clip_vert2 = clip_polygon.vertices[(clip_i > 0) ? (clip_i - 1) : clip_n - 1];
            Vector cur_vert = subject_polygon.vertices[subj_i];
            Vector prv_vert = subject_polygon.vertices[(subj_i > 0) ? (subj_i - 1) : subj_n - 1];

            Vector intersection = intersect(prv_vert, cur_vert, clip_vert1, clip_vert2);

            if (is_inside(cur_vert, clip_vert1, clip_vert2, clip_center)) {
                if (!is_inside(prv_vert, clip_vert1, clip_vert2, clip_center)) {
                    out_polygon.vertices.push_back(intersection);
                }
                out_polygon.vertices.push_back(cur_vert);
            }

            else if (is_inside(prv_vert, clip_vert1, clip_vert2, clip_center)) {
                out_polygon.vertices.push_back(intersection);
            }
        }

        subject_polygon = out_polygon;

    }

    return out_polygon;

}

Polygon clip_polygon(Polygon subject_polygon, Vector &clip_vert1, Vector &clip_vert2, Vector &site) {

    Polygon out_polygon = Polygon();
    int subj_n = subject_polygon.vertices.size();

    for (int subj_i = 0; subj_i < subj_n; ++subj_i) {

        Vector cur_vert = subject_polygon.vertices[subj_i];
        Vector prv_vert = subject_polygon.vertices[(subj_i > 0) ? (subj_i - 1) : subj_n - 1];

        Vector intersection = intersect(prv_vert, cur_vert, clip_vert1, clip_vert2);

        if (is_inside(cur_vert, clip_vert1, clip_vert2, site)) {
            if (!is_inside(prv_vert, clip_vert1, clip_vert2, site)) {
                out_polygon.vertices.push_back(intersection);
            }
            out_polygon.vertices.push_back(cur_vert);
        }

        else if (is_inside(prv_vert, clip_vert1, clip_vert2, site)) {
            out_polygon.vertices.push_back(intersection);
        }
    }

    return out_polygon;

}


vector<Polygon> voronoi(vector<Vector> points, Polygon site_start_polygon, vector<double> weights) {

    vector<Polygon> voronoi_diagram;

    // #pragma omp parallel
    for(int i = 0; i < points.size(); ++i) {
        
        Vector Pi = points[i];
        Polygon site_polygon = site_start_polygon;

        for(int j = 0; j < points.size(); ++j) {
            
            Vector Pj = points[j];

            if (norm(Pi - Pj) == 0)
                continue;

        Vector M;
        if (!weights.size())
            M = (Pj + Pi) / 2;
        else
            M = (Pj + Pi) / 2 + ((weights[i] - weights[j]) / (2 * pow(norm(Pi - Pj), 2))) * (Pj - Pi);
    
        Vector bissector_point = M + Vector(Pi[1] - Pj[1], Pj[0] - Pi[0], 0);

        site_polygon = clip_polygon(site_polygon, M, bissector_point, Pi);     

        }
        
        voronoi_diagram.push_back(site_polygon);

    }

    return voronoi_diagram;

}
