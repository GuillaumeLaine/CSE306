#include <vector>
#include <cmath>
#include <math.h>
#include <random>
#include <algorithm>

#include "objects.h"
#include "mesh.cpp"

using namespace std;

static default_random_engine engine(10);
static uniform_real_distribution<double> uniform(0, 1);

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
    return scale(1/norm(a), a);
}

Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}


// Ray definition
Ray::Ray(Vector origin, Vector direction) {
    O = origin;
    u = direction;
};

// Sphere definition
Sphere::Sphere(Vector center, double radius, Vector color, bool reflects, bool refracts) {
    C = center;
    R = radius;
    albedo = color;
    this->reflects = reflects;
    this->refracts = refracts;
};

Intersection Sphere::intersect(const Ray& r) {
    Intersection i;
    Vector v = r.O - C;
    double d = dot(r.u, v);
    
    double delta = d * d - (dot(v, v) - (R * R));
    
    if (delta >= 0) {
    
        double sq_delta = sqrt(delta);
        double minus_d = dot(r.u, -1 * v);
        double t2 = minus_d + sq_delta;
    
        if (t2 >= 0) {
    
            double t1 = minus_d - sq_delta;
    
            if (t1 >= 0) {
                i.flag = true;
                i.P = r.O + t1 * r.u;
                i.dist_t = t1;
            }
    
            else {
                i.flag = true;
                i.P = r.O + t2 * r.u;
                i.dist_t = t2;
            }

            i.N = unit(i.P - C);
            i.albedo = albedo;
        }
    }

    return i;
}

// Scene definition
Scene::Scene(vector<Geometry*> spheres) {
    this->s = spheres;
};

Intersection Scene::intersect(const Ray& r) {
    
    Intersection intersect_point;
    double min_dist = 100000000;

    for (unsigned int i=0; i<s.size(); i++) {

        Geometry* sphere = s[i];
        Intersection intersection = sphere->intersect(r);

        if (intersection.flag) {

            double dist_sq = dot(r.O - intersection.P, r.O - intersection.P); // skip sqrt computation for rapidity

            if (dist_sq <= min_dist) {

                min_dist = dist_sq;
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

    if (i.flag) {
    
        Geometry* int_sphere = s[i.sphere_id];

        // Reflective surface
        if (int_sphere->reflects) {
            
            Ray reflected_ray(i.P + scale(0.01, i.N), r.u - (scale(2 * dot(r.u, i.N), i.N)));
            return getColor(reflected_ray, S, ray_depth - 1);
        
        }

        // Transparent surface
        if (int_sphere->refracts) {

            const double glass_n = 1.5168;

            bool is_entering = dot(r.u, i.N) < 0;  // check ray is entering ball
            double n1 = is_entering ? 1 : glass_n;
            double n2 = is_entering ? glass_n : 1;
            Vector N = is_entering ? i.N : -1 * i.N;

            double k0 = ((n1 - n2) * (n1 - n2)) / ((n1 + n2) * (n1 + n2));
            double R = k0 + (1 - k0) * pow(1 - abs(dot(N, r.u)), 5);

            double rand_float = (double) (rand() % 100) / 100;

            // reflect with probability, to model reflecting portion of ray
            if (rand_float < R) {

                Ray reflected_ray(i.P + 0.01*N, r.u - (scale(2 * dot(r.u, N), N)));
                return getColor(reflected_ray, S, ray_depth - 1);

            }

            // otherwise refract "the rest of the ray"
            Vector wT = (n1 / n2) * (r.u - dot(r.u, N) * N);
            Vector wN = -1 * N * (sqrt(1 - (n1/n2) * (n1/n2) * (1 - pow(dot(r.u, N), 2))));
            Vector w = wT + wN;
            Ray refracted_ray(i.P - 0.01*N, w); // in refraction, move point slightly to the other side

            return getColor(refracted_ray, S, ray_depth - 1);

        }

        // Diffuse surface
        else {
            
            Vector direct_contribution = intensity(*this, i, S);
        
            Ray random_ray = Ray(i.P + 0.01*i.N, random_cos(i.N));
            Vector indirect_contribution = i.albedo * getColor(random_ray, S, ray_depth - 1);

            return (direct_contribution + indirect_contribution) / 2;

        }
    }

    // No ray intersection
    return Vector();
}


// Utils
Vector intensity(Scene& scene, Intersection& i, Vector& S) {
    double I = 200000000;
    Vector v1 = S - i.P;
    double d_sq = dot(v1, v1);
    Vector wi = v1 / sqrt(d_sq);
    
    double vp = 0;

    //Check visibility
    Ray r(i.P + 0.01*i.N, wi);
    Intersection light_i = scene.intersect(r);

    if (!light_i.flag) {
        vp = 1; //visible    
    }

    else {
        double v2 = light_i.dist_t;
        if (v2 * v2 > d_sq) {
            vp = 1;
        }
    }

    return i.albedo * (I * vp * max(dot(i.N, wi), 0.) / (4 * M_PI * M_PI * d_sq));
};

Vector random_cos(const Vector &N) {

    double r1 = uniform(engine);
    double r2 = uniform(engine);

    double x = cos(2 * M_PI * r1) * sqrt(1 - r2);
    double y = sin(2 * M_PI * r1) * sqrt(1 - r2);
    double z = sqrt(r2);

    vector<double> v = {abs(N[0]), abs(N[1]), abs(N[2])};

    int minElemInd = min_element(v.begin(), v.end()) - v.begin();

    Vector T1(0, 0, 0);
    T1[minElemInd] = 0;
    T1[(minElemInd + 1) % 3] = - N[(minElemInd - 1) % 3];
    T1[(minElemInd - 1) % 3] = N[(minElemInd + 1) % 3];
    T1 = unit(T1);

    Vector T2 = cross(N, T1);

    return x * T1 + y * T2 + z * N;
}

Vector boxMuller(double stdev) {

    double r1 = uniform(engine);
    double r2 = uniform(engine);
    double x = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2) * stdev;
    double y = sqrt(-2 * log(r1)) * sin(2 * M_PI * r2) * stdev;

    return Vector(x, y, 0.);

}
