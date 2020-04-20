#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image_write.h"
#include "objects.h"

using namespace std;


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

            Vector color = scene.getColor(r, S);
            img[(i * W + j)*3 + 0] = max(0., min(255., pow(color[0], gamma)));
            img[(i * W + j)*3 + 1] = max(0., min(255., pow(color[1], gamma)));
            img[(i * W + j)*3 + 2] = max(0., min(255., pow(color[2], gamma)));

        }
    }

    // Write to image
    stbi_write_png("image.png", W, H, 3, &img[0], 0);

    return 0;
}
