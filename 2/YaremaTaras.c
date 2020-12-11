#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX_ITERATIONS 1000
#define TOLERANCE 1e-8

/*
 * Returns the value f(x, y) for (x, y) in R^2 
 */
double f(double x, double y) {
    return - 8 - 36 * x - 40 * pow(x, 2) + 18 * pow(x, 3) + 9 * pow(x, 4) - 36 * y + 18 * pow(x, 2) * y
        + 20 * pow(y, 2) + 32 * x * pow(y, 2) - 21 * pow(x, 2) * pow(y, 2) - 2 * pow(x, 3) * pow(y, 2) 
        + 5 * pow(x, 4) * pow(y, 2) + 5 * pow(y, 3) - 2 * pow(x, 2) * pow(y, 3) - 30 * pow(y, 4) 
        - 2 * x * pow(y, 4) + 30 * pow(x, 2) * pow(y, 4) - 2 * pow(y, 5) + 4 * pow(y, 6);
}

/*
 * Computes the partial derivative of f respect to x
 * and evaluates in (x, y)
 */
double partial_x(double x, double y) {
    return -36 - 80 * x + 54 * pow(x, 2) + 36 * pow(x, 3) + 36 * x * y + 32 * pow(y, 2) - 42 * x * pow(y, 2)
        - 6 * pow(x, 2) * pow(y, 2) + 20 * pow(x, 3) * pow(y, 2) - 4 * x * pow(y, 3) - 2 * pow(y, 4) + 60 * x * pow(y, 4);
}

/*
 * Computes the partial derivative of f respect to y
 * and evaluates in (x, y)
 */
double partial_y(double x, double y) {
    return - 36 + 18 * pow(x, 2) +40 * y + 64 * x * y - 42 * pow(x, 2) * y - 4 * pow(x, 3) * y + 10 * pow(x, 4) * y 
        + 15 * pow(y, 2) - 6 * pow(x, 2) * pow(y, 2) - 120 * pow(y, 3) - 8 * x * pow(y, 3) + 120 * pow(x, 2) * pow(y, 3) 
        - 10 * pow(y, 4) + 24 * pow(y, 5);
}

/*
 * Computes the gradient of f and evaluates in (x, y)
 */
void gradient(double x, double y, double *r) {
    r[0] = partial_x(x, y);
    r[1] = partial_y(x, y);
}

/*
 * Newton method to find solutions with the y component zero
 */
int newton_zero_x(double *x0) {
    int its = 0;
    double x;

    while (its < MAX_ITERATIONS) {
        x = *x0 - f(*x0, 0) / partial_x(*x0, 0);
        if (fabs(f(*x0, 0)) < TOLERANCE) {
            break;
        }

        *x0 = x;
        its++;
    }

    return its;
}

/*
 * Newton method to find solutions with the x component zero
 */
int newton_zero_y(double *y0) {
    int its = 0;
    double y;

    while (its < MAX_ITERATIONS) {
        y = *y0 - f(0, *y0) / partial_y(0, *y0);
        if (fabs(f(0, *y0)) < TOLERANCE) {
            break;
        }

        *y0 = y;
        its++;
    }

    return its;
}

int main() {
    int its_x, its_y;
    double x, y;

    x = 0;
    y = 0;

    its_x = newton_zero_x(&x);
    its_y = newton_zero_y(&y);

    printf("f(%f, 0) -> (%f, 0) in %d its\n", 1.0, x, its_x);
    printf("f(0, %f) -> (0, %f) in %d its\n", 1.0, y, its_y);

    return 0;
}
