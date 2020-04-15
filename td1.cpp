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

    void print_coords() {
        cout << coords[0] << ' ' << coords[1] << ' ' << coords[2] << '\n';
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

double norm(const Vector& a) {
    return sqrt(pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2));
}

Vector unit(const Vector& a) {
    return scale(1/norm(a), a);
}

double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
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
        
        double delta = pow(d, 2) - (pow(norm(r.O - C), 2) - (R * R));
        if (delta >= 0) {
            double sq_delta = sqrt(delta);
            double t1 = d - sq_delta;
            double t2 = d + sq_delta;
            if (t2 >= 0) {
                if (t1 >= 0) {
                    i.flag = true;
                    i.P = r.O + scale(t1, r.u);
                }
                else {
                    i.flag = true;
                    i.P = r.O + scale(t2, r.u);
                }
                Vector n = i.P - C;
                i.N = unit(n);
            }
            i.albedo = albedo;
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
                Vector cam_pos = Vector(0, 0, 55);
                double dist = norm(cam_pos - intersection.P);
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
    Ray r(i.P + scale(0.01*norm(wi), wi), wi);
    Intersection light_i = scene.intersect(r);

    if (!light_i.flag) {
        vp = 1; //visible    
    }

    else {
        if (norm(light_i.P - S) >= d) {
            vp = 1;
        }
    }

    return (I * vp * max(-dot(i.N, wi), 0.) / (4 * M_PI * M_PI * d * d));
};


int main() {

    // Define scene
    Sphere s(Vector(0, 0, 0), 10, Vector(1, 0, 1));
    /*
    Sphere s1(Vector(0, 1000, 0), 940, Vector(1, 0, 0)); 
    Sphere s2(Vector(0, 0, -1000), 940, Vector(0, 1, 0));
    Sphere s3(Vector(0, -1000, 0), 990, Vector(0, 0, 1));
    Sphere s4(Vector(0, 0, 1000), 940, Vector(1, 0, 1));
    */

    vector<Sphere> scene_vector;
    scene_vector.push_back(s);
    /*
    scene_vector.push_back(s1);
    scene_vector.push_back(s2);
    scene_vector.push_back(s3);
    scene_vector.push_back(s4);
    */

    Scene scene(scene_vector);
    Vector S(-10, 20, 40);

    // Viewing window
    int H = 400;
    int W = 400;
    Vector cam(0, 0, 55);
    double fov = M_PI / 3;
    double z = cam[2] - W/(2*tan(fov/2));

    // Ray trace
    vector<unsigned char> img;
    for (int i = 0; i < W; i++) {
        for (int j = 0; j < H; j++) {
            double x = j;
            double y = H - i - 1;
            Vector pixel(
                cam[0] + x + 0.5 - W/2,
                cam[1] + y + 0.5 - H/2,
                z
            );
            Vector ray_dir = pixel - cam;
            ray_dir = unit(ray_dir);
            Ray r(pixel, ray_dir);

            Intersection intersect = scene.intersect(r);
            if (intersect.flag) {
                Vector albedo = intersect.albedo;
                double L = intensity(scene, intersect, S);
                img.push_back(albedo[0] * L * 255);
                img.push_back(albedo[1] * L * 255);
                img.push_back(albedo[2] * L * 255);
            }
            else {
                img.push_back(0);
                img.push_back(0);
                img.push_back(0);
            }

        }
    }
    stbi_write_png("image.png", W, H, 3, &img[0], 0);

    return 0;
}
