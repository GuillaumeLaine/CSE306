#include <random>
#include <math.h>
#include <iostream>
#define PI 3.14159265

//Monte Carlo integration example
//output : 24.7429 compared to real value of 24.3367

static std::default_random_engine engine(10); 
static std::uniform_real_distribution<double> uniform(0, 1);

double p(double x, double y, double z) {
    return pow((1 / sqrt(2*PI)), 3) * exp(-(x*x + y*y + z*z) / 2);
}

double f(double x, double y, double z) {
    return cos(x*y*z);
}

int main() {
    double sum = 0;
    for(int i = 0; i < 10000; i++) {
        double r1 = uniform(engine);
        double r2 = uniform(engine);
        double r3 = uniform(engine);
        double r4 = uniform(engine);
        double x = sqrt(-2 * log(r1)) * cos(2 * PI * r2);
        double y = sqrt(-2 * log(r1)) * sin(2 * PI * r2);
        double z = sqrt(-2 * log(r3)) * cos(2 * PI * r4);
        sum += f(x, y, z) / p(x, y, z);
    }
    std::cout << sum/10000;
    return 0;
}