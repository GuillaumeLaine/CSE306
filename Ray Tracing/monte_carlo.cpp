#include <random>
#include <math.h>
#include <iostream>

static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0, 1);

double f(double x, double y, double z) {
    return cos(x * y * z);
}

double gauss(double x, double y, double z) {
    return pow(1 / sqrt(2 * M_PI), 3) * exp(- (x * x + y * y + z * z) / 2 );
}

int main() {

    double S = 0;
    int it = 10000;

    for (int i = 0; i < it; i++) {
        
        double r1 = uniform(engine);
        double r2 = uniform(engine);
        double r3 = uniform(engine);
        double r4 = uniform(engine);

        double x = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2);
        double y = sqrt(-2 * log(r1)) * sin(2 * M_PI * r2);
        double z = sqrt(-2 * log(r3)) * cos(2 * M_PI * r4); // idea for z formula from https://stats.stackexchange.com/questions/51321/how-to-use-box-muller-transform-to-generate-n-dimensional-normal-random-variable
        
        S += f(x, y, z) / gauss(x, y ,z);
    }

    std::cout << "Estimation is: " << S / it;

    return 0;

}