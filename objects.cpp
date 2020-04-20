#include <vector>
#include <cmath>
#include <math.h>

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
    return  Vector(a * v[0], a * v[1], a * v[2]);
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
    return scale(1/norm(a), a);
}


// Ray definition
Ray::Ray(Vector origin, Vector direction) {
    O = origin;
    u = direction;
};

// Sphere definition
Sphere::Sphere(Vector center, double radius, Vector color, bool reflects) {
    C = center;
    R = radius;
    albedo = color;
    this->reflects = reflects; 
};

Intersection Sphere::intersect(const Ray& r) {
    Intersection i;
    double d = dot(r.u, r.O - C);
    
    double delta = d * d - (pow(norm(r.O - C), 2) - (R * R));
    if (delta >= 0) {
        double sq_delta = sqrt(delta);
        double minus_d = dot(r.u, C - r.O);
        double t1 = minus_d - sq_delta;
        double t2 = minus_d + sq_delta;
        if (t2 >= 0) {
            if (t1 >= 0) {
                i.flag = true;
                i.P = r.O + scale(t1, r.u);
            }
            else {
                i.flag = true;
                i.P = r.O + scale(t2, r.u);
            }

            i.N = unit(i.P - C);
            i.albedo = albedo;
        }
    }

    return i;
}

// Scene definition
Scene::Scene(vector<Sphere> spheres) {
    this->s = spheres;
};

Intersection Scene::intersect(const Ray& r) {
    
    Intersection intersect_point;
    double min_dist = 1000000;

    for (int i=0; i<s.size(); i++) {
        Sphere sphere = s[i];
        Intersection intersection = sphere.intersect(r);
        if (intersection.flag) {
            double dist = norm(r.O - intersection.P);
            if (dist <= min_dist) {
                min_dist = dist;
                intersect_point = intersection;
                intersect_point.sphere_id = i;
            }
        }
    }

    return intersect_point;
}

Vector Scene::getColor(const Ray& r, Vector& S, int ray_depth) {

    if (ray_depth < 0) {
        return Vector(0., 0., 0.);
    }

    Intersection i = this->intersect(r);
    int sphere_id = i.sphere_id;

    if (i.flag) {

        if (s[sphere_id].reflects) {
            
            Ray reflected_ray(i.P + scale(0.01, i.N), r.u - (scale(2 * dot(r.u, i.N), i.N)));
            return getColor(reflected_ray, S, ray_depth - 1);
        
        }

        else {
        
            return intensity(*this, i, S);
        
        }
    }
}


// Utils
Vector intensity(Scene& scene, Intersection& i, Vector& S) {
    double I = 100000000;
    double d = norm(S - i.P);
    Vector wi = unit(S - i.P);
    int vp = 0;

    //Check visibility
    Ray r(i.P + scale(0.01, i.N), wi);
    Intersection light_i = scene.intersect(r);

    if (!light_i.flag) {
        vp = 1; //visible    
    }

    else {
        if (norm(light_i.P - i.P) > d) {
            vp = 1;
        }
    }

    return i.albedo * (I * vp * max(dot(i.N, wi), 0.) / (4 * M_PI * M_PI * d * d)) * 255;
};
