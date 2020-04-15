#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
using namespace std;

class Vector {

public:
    explicit Vector(double x = 0., double y = 0., double z = 0.) {
        coords[0] = x;
        coords[1] = y;
        coords[2] = z;
    };

    Vector& operator+=(const Vector& b) {
        coords[0] += b[0];
        coords[1] += b[1];
        coords[2] += b[2];
        return *this;
    }

    const double& operator[](int i) const {return coords[i]; }
    double& operator[](int i) { return coords[i]; }

private: 
    double coords[3];

};

Vector operator+(const Vector& a, const Vector &b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector operator-(const Vector& a, const Vector &b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
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

struct Intersection {

    Vector P;
    Vector N;
    Vector albedo; 
    bool flag = false;
};

class Ray {

public:
    explicit Ray(Vector origin, Vector direction) {
        O = origin;
        u = direction;
    };

    Vector O;
    Vector u;

};

class Sphere {

public:
    explicit Sphere(Vector center, double radius, Vector color) {
        C = center;
        R = radius;
        albedo = color;
    };

    Intersection intersect(Ray& r) {
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

    double R;
    Vector C;
    Vector albedo; 

};

class Scene {

public:
    Scene(vector<Sphere> spheres) {
        this->s = spheres;
    }

    Intersection intersect(Ray& r) {
        
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
                }
            }
        }

        return intersect_point;
    }

    vector<Sphere> s;

};

struct TriangleIndices {

    int vtxindices[3];
    int normalindices[3];
    int uvindices[3];

};

double intensity(Scene& scene, Intersection& i, Vector& S) {
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

    return (I * vp * max(dot(i.N, wi), 0.) / (4 * M_PI * M_PI * d * d));
};


int main() {

    // Define scene
    Sphere s(Vector(0, 0, 0), 10, Vector(1, 1, 1));
    Sphere s1(Vector(0, 1000, 0), 940, Vector(1, 0, 0)); 
    Sphere s2(Vector(0, 0, -1000), 940, Vector(0, 1, 0));
    Sphere s3(Vector(0, -1000, 0), 990, Vector(0, 0, 1));
    Sphere s4(Vector(0, 0, 1000), 940, Vector(1, 0, 1));
    Sphere s5(Vector(-1000, 0, 0), 940, Vector(0, 1, 1));
    Sphere s6(Vector(1000, 0, 0), 940, Vector(1, 1, 0));

    vector<Sphere> scene_vector;
    scene_vector.push_back(s);
    scene_vector.push_back(s1);
    scene_vector.push_back(s2);
    scene_vector.push_back(s3);
    scene_vector.push_back(s4);
    scene_vector.push_back(s5);
    scene_vector.push_back(s6);

    Scene scene(scene_vector);
    Vector S(-10, 20, 40);

    // Viewing window
    int H = 512;
    int W = 512;
    Vector cam(0, 0, 55);
    double fov = M_PI / 3;
    double z = cam[2] - W/(2*tan(fov/2));
    double gamma = 1. / 2.2;

    // Ray trace
    vector<unsigned char> img(W*H*3);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            double x = j;
            double y = H - i - 1;
            Vector pixel(
                cam[0] + x + 0.5 - W/2,
                cam[1] + y + 0.5 - H/2,
                z
            );
            Vector ray_dir = unit(pixel - cam);
            Ray r(cam, ray_dir);

            Intersection intersect = scene.intersect(r);
            if (intersect.flag) {
                Vector albedo = intersect.albedo;
                double L = intensity(scene, intersect, S);
                double throuput = L * 255;
                img[(i * W + j)*3 + 0] = max(0., min(255., pow(albedo[0] * throuput, gamma)));
                img[(i * W + j)*3 + 1] = max(0., min(255., pow(albedo[1] * throuput, gamma)));
                img[(i * W + j)*3 + 2] = max(0., min(255., pow(albedo[2] * throuput, gamma)));
            }
            else {
                img[(i * W + j)*3 + 0] = 0;
                img[(i * W + j)*3 + 1] = 0;
                img[(i * W + j)*3 + 2] = 0;
            }

        }
    }
    stbi_write_png("image.png", W, H, 3, &img[0], 0);

    return 0;
}
