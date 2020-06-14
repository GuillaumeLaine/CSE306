#include "objects.h"
#include "utils.cpp"

int main() {

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

    save_svg({subject}, "polygons/subject.svg");
    save_svg({clipper}, "polygons/clipper.svg");

    Polygon clipped = clip_polygon(subject, clipper);

    save_svg({clipped}, "polygons/clipped.svg");

    return 0;

}