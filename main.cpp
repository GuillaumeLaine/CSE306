#include "objects.h"
#include "utils.cpp"
#include "time.h"

void test_poly_clipping() {

    // clipping polygon test
    Polygon subject = Polygon({
        Vector(0.2, 0.6, 0),
        Vector(0.8, 0.6, 0), 
        Vector(0.7, 0.5, 0), 
        Vector(0.8, 0.4, 0), 
        Vector(0.2, 0.4, 0), 
        Vector(0.3, 0.5, 0)
        });

    Polygon clipper = Polygon({
        Vector(0.9, 0.9, 0), 
        Vector(0.9, 0.1, 0), 
        Vector(0.2, 0.1, 0)});

    save_svg({subject}, "results/subject.svg");
    save_svg({clipper}, "results/clipper.svg");

    Polygon clipped = clip_polygon(subject, clipper);

    save_svg({clipped}, "results/clipped.svg");

}

void test_voronoi() {

    std::vector<Vector> points = {
        Vector(0.2, 0.6, 0),
        Vector(0.8, 0.6, 0), 
        Vector(0.7, 0.5, 0), 
        Vector(0.8, 0.4, 0), 
        Vector(0.2, 0.4, 0), 
        Vector(0.3, 0.5, 0)
        };

    Polygon space({
        Vector(0, 0, 0),
        Vector(0, 1, 0), 
        Vector(1, 1, 0), 
        Vector(1, 0, 0)});

    std::vector<Polygon> diagram = voronoi(points, space);
    save_svg_with_points(diagram, points, "results/voronoi.svg");

    std::vector<Vector> more_points = generate_points(10000);
    save_svg_with_points(voronoi(more_points, space), more_points, "results/voronoi_more.svg");

}
int main() {

    clock_t tStart = clock();

    test_poly_clipping();
    test_voronoi();

    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    return 0;

}