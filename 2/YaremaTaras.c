#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct point {
    double x;
    double y;
} Point;


double f(double x, double y) {
    double z;
    
    z = -8 - 36 * x - 40 * pow(x, 2) + 18 * pow(x, 3) + 9 * pow(x, 4) - 36 * y + 18 * pow(x, 2) * y;
    z += 20 * pow(y, 2) + 32 * x * pow(y, 2) - 21 * pow(x, 2) * pow(y, 2) - 2 * pow(x, 3) * pow(y, 2) + 5 * pow(x, 4) * pow(y, 2);
    z += 5 * pow(y, 3) - 2 * pow(x, 2) * pow(y, 3) - 30 * pow(y, 4) - 2 * x * pow(y, 4) + 30 * pow(x, 2) * pow(y, 4) - 2 * pow(y, 5) + 4 * pow(y, 6);

    return z;
}

Point gf(double x, double y) {
    Point gradient = {0, 0};

    gradient.x = -36 - 80 * x + 54 * pow(x, 2) + 36 * pow(x, 3) + 36 * x * y + 32 * pow(y, 2) - 42 * x * pow(y, 2);
    gradient.x += -6 * pow(x, 2) * pow(y, 2) + 20 * pow(x, 3) * pow(y, 2) - 4 * x * pow(y, 3) - 2 * pow(y, 4) + 60 * x * pow(y, 4);

    gradient.y = -36 + 18 * pow(x, 2) +40 * y + 64 * x * y - 42 * pow(x, 2) * y - 4 * pow(x, 3) * y + 10 * pow(x, 4) * y + 15 * pow(y, 2) - 6 * pow(x, 2) * pow(y, 2);
    gradient.y += -120 * pow(y, 3) - 8 * x * pow(y, 3) + 120 * pow(x, 2) * pow(y, 3) - 10 * pow(y, 4) + 24 * pow(y, 5);

    return gradient;
}

int main() {
    int i, j;
    for (i = 0; i <= 100; i++) {
        for (j = 0; j <= 100; j++) {
            printf("f(%f, %f) = %f\n", (double)i, (double)j, f((double)i, (double)j));
            Point grad = gf((double)i, (double)j);
            printf("∇f(%f, %f) = {%f, %f}\n", (double)i, (double)j, grad.x, grad.y);
        }
    }
    return 0;
}
