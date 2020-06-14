#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include <time.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image_write.h"
#include "mesh.cpp"
#include "objects.h"

using namespace std;


int main() {

    clock_t tStart = clock();

    // Define scene
    Sphere* white_ball  = new Sphere(Vector(0, 0, 0), 10, Vector(1, 1, 1));
    Sphere* ceiling     = new Sphere(Vector(0, 1000, 0), 940, Vector(1, 0, 0)); 
    Sphere* front_wall  = new Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0));
    Sphere* floor       = new Sphere(Vector(0, -1000, 0), 990, Vector(0, 0, 1));
    Sphere* back_wall   = new Sphere(Vector(0, 0, 1000), 940, Vector(1, 0, 1));
    Sphere* left_wall   = new Sphere(Vector(-1000, 0, 0), 940, Vector(0, 1, 1));
    Sphere* right_wall  = new Sphere(Vector(1000, 0, 0), 940, Vector(1, 1, 0));
    Sphere* mirror_ball = new Sphere(Vector(-20, 0, 0), 10, Vector(1, 1, 1));
    Sphere* glass_ball  = new Sphere(Vector(20, 0, 0), 10, Vector(1, 1, 1));
    mirror_ball->reflects = true;
    glass_ball->refracts = true;

    TriangleMesh* cat = new TriangleMesh();
    cat->readOBJ("./cadnav/Models_F0202A090/cat.obj");
    for (unsigned int i=0; i<cat->vertices.size(); ++i) {
        cat->vertices[i] = cat->vertices[i] * 0.6 + Vector(0, -10, 0);
    }

    Scene scene({
        cat,
        // white_ball, mirror_ball, glass_ball,
        front_wall, back_wall, left_wall, right_wall, ceiling, floor
        });
    Vector S(-10, 20, 40);

    // Viewing window
    const int H = 512;
    const int W = 512;
    Vector cam(0, 0, 55);
    double fov = M_PI / 3;
    double z = cam[2] - W/(2*tan(fov/2));
    double gamma = 1. / 2.2;
    const int K = 15 ;

    // Ray trace
    vector<unsigned char> img(W*H*3);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            
            double x = j;
            double y = H - i - 1;

            Vector pixel(cam[0] + x + 0.5 - W/2, cam[1] + y + 0.5 - H/2, z);

            // Average color over each pixel ray
            Vector color = Vector();
    
            for (int k = 0; k < K; k++) {

                Vector pixel_aa = pixel + boxMuller(); // offset pixel positioning for anti-aliasing                 
                Vector ray_dir = unit(pixel_aa - cam);
                Ray r(cam, ray_dir);
                color += scene.getColor(r, S);

            }
            color = color * 255 / K;

            img[(i * W + j)*3 + 0] = max(0., min(255., pow(color[0], gamma)));
            img[(i * W + j)*3 + 1] = max(0., min(255., pow(color[1], gamma)));
            img[(i * W + j)*3 + 2] = max(0., min(255., pow(color[2], gamma)));

        }
    }

    // Write to image
    stbi_write_png("image.png", W, H, 3, &img[0], 0);

    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}

