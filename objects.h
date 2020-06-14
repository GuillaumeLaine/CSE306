#pragma once

#include <vector>
#include <cmath>
#include <math.h>

using namespace std;

class Vector {

public:
    explicit Vector(double x = 0., double y = 0., double z = 0.);
    Vector& operator+=(const Vector& b);
    const double& operator[](int i) const;
    double& operator[](int i);

private: 
    double coords[3];
};

class Polygon {

public:
    vector<Vector> vertices;

    Polygon() {};

    Polygon(vector<Vector> vertices) {
        this->vertices = vertices;
    }

    Vector center();

};

Vector operator+(const Vector& a, const Vector &b);

Vector operator-(const Vector& a, const Vector &b);

Vector operator*(double a, const Vector& v);

Vector operator*(const Vector& v, double a);

Vector operator/(const Vector& v, double a);

Vector scale(double a, const Vector& v);

double dot(const Vector& a, const Vector& b);

double norm(const Vector& a);

Vector unit(const Vector& a);

Vector cross(const Vector& a, const Vector& b);

bool is_inside(Vector const &P, Vector const &u, Vector const &v, Vector const &poly_center);

// Check if [AB] intersects (uv)
Vector intersect(Vector const &A, Vector const &B, Vector const &u, Vector const &v);

Polygon clip_polygon(Polygon subject_polygon, Polygon &clip_polygon);

Polygon clip_polygon(Polygon subject_polygon, Vector &clip_vert1, Vector &clip_vert2, Vector &site);

vector<Polygon> voronoi(vector<Vector> points, Polygon site_start_polygon);