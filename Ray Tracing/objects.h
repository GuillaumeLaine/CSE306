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

struct Intersection {

    Vector P;
    Vector N;
    Vector albedo;
    bool flag = false;
    int sphere_id;
    double dist_t;
};


class Ray {

public:
    Vector O;
    Vector u;

    explicit Ray(Vector origin, Vector direction);
};

class Geometry {

public:
    
    Vector albedo = Vector(1, 1, 1);
    bool reflects;
    bool refracts;

    virtual Intersection intersect (const Ray& r) = 0;

};

class Sphere : public Geometry{

public:

    double R;
    Vector C;

    explicit Sphere(Vector center, double radius, Vector color, bool reflects = false, bool refracts = false);
    Intersection intersect(const Ray& r);
};

class Scene {

public:

    vector<Geometry*> s;
    
    Scene(vector<Geometry*> spheres);
    Intersection intersect(const Ray& r);
    Vector getColor(const Ray& r, Vector& S,  int ray_depth = 5);

};

Vector intensity(Scene& scene, Intersection& i, Vector& S);

Vector random_cos(const Vector &N);

Vector boxMuller(double stdev=1.);