#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define TOLERANCE 1e-12
#define H_STEP 1e-3
#define MAX_ITERATIONS 1000
#define PLOT_STEPS 100000

/*
 * Computes the dot product of the x and y
 */
double dot(double *x, double *y) { return x[0] * y[0] + x[1] * y[1]; }

/*
 * Returns the value f(x, y) for (x, y) in R^2
 */
double f(double x, double y) {
    return -8 - 36 * x - 40 * pow(x, 2) + 18 * pow(x, 3) + 9 * pow(x, 4) -
           36 * y + 18 * pow(x, 2) * y + 20 * pow(y, 2) + 32 * x * pow(y, 2) -
           21 * pow(x, 2) * pow(y, 2) - 2 * pow(x, 3) * pow(y, 2) +
           5 * pow(x, 4) * pow(y, 2) + 5 * pow(y, 3) -
           2 * pow(x, 2) * pow(y, 3) - 30 * pow(y, 4) - 2 * x * pow(y, 4) +
           30 * pow(x, 2) * pow(y, 4) - 2 * pow(y, 5) + 4 * pow(y, 6);
}

/*
 * Computes the partial derivative of f respect to x
 * and evaluates in (x, y)
 */
double partial_x(double x, double y) {
    return -36 - 80 * x + 54 * pow(x, 2) + 36 * pow(x, 3) + 36 * x * y +
           32 * pow(y, 2) - 42 * x * pow(y, 2) - 6 * pow(x, 2) * pow(y, 2) +
           20 * pow(x, 3) * pow(y, 2) - 4 * x * pow(y, 3) - 2 * pow(y, 4) +
           60 * x * pow(y, 4);
}

/*
 * Computes the partial derivative of f respect to y
 * and evaluates in (x, y)
 */
double partial_y(double x, double y) {
    return -36 + 18 * pow(x, 2) + 40 * y + 64 * x * y - 42 * pow(x, 2) * y -
           4 * pow(x, 3) * y + 10 * pow(x, 4) * y + 15 * pow(y, 2) -
           6 * pow(x, 2) * pow(y, 2) - 120 * pow(y, 3) - 8 * x * pow(y, 3) +
           120 * pow(x, 2) * pow(y, 3) - 10 * pow(y, 4) + 24 * pow(y, 5);
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

/*
 * Computes the tangent vector to the curve
 * in the coordinates p with norm 1
 *
 * Returns an integer as we may encounter divisions by zero
 * (non-zero return in that case)
 */
int tangent(double *p, double *vec) {
    double norm;

    vec[0] = -partial_y(p[0], p[1]);
    vec[1] = partial_x(p[0], p[1]);

    norm = sqrt(pow(vec[0], 2) + pow(vec[1], 2));

    /* Check that we do not divide by zero */
    if (norm == 0) {
        return 1;
    }

    /* Normalize to have norm 1 */
    vec[0] = vec[0] / norm;
    vec[1] = vec[1] / norm;

    return 0;
}

/*
 * Uses the Newton method to find a zero of f at distance h from a
 * first zero (x0, y0)
 */
int newton_iterator(double x0, double y0, double x1, double y1, double h,
                    double *result) {
    double *H, *partial_f;
    double detH, x2, y2, incr_x, incr_y;
    int n = 0, exit_code = 0;

    H = calloc(2, sizeof(double));
    if (H == NULL) return 2;

    partial_f = calloc(2, sizeof(double));
    if (partial_f == NULL) return 2;

    /* Loop for MAX_ITERATIONS max times */
    while (n < MAX_ITERATIONS) {
        /*
         * We will stop when both of the following conditions are met,
         * defined by the system in problem 3:
         *  1. The absolute value of f(x_k, y_k) is sufficiently near (0, 0)
         *  2. The value (x_k - x_0)^2 + (y_k - y_0)^2 = h^2 is sufficiently low
         */
        if (fabs(f(x1, y1)) < TOLERANCE &&
            fabs(pow(x1 - x0, 2) + pow(y1 - y0, 2) - pow(h, 2)) < TOLERANCE) {
            break;
        }

        /* Compute H(x_k) = (f(x_k), g(x_k)) */
        H[0] = f(x1, y1);
        H[1] = pow(x1 - x0, 2) + pow(y1 - y0, 2) - pow(h, 2);

        /* Compute the determinant of DH */
        gradient(x1, y1, partial_f);
        detH = 2 * (y1 - y0) * partial_f[0] - 2 * (x1 - x0) * partial_f[1];

        /* We got determinant zero, so we exit with error */
        if (detH == 0) {
            exit_code = 1;
            goto out;
        }

        /* Compute the product DH(x_k)^{-1} * H(x_k)^T */
        incr_x = (2 * (y1 - y0) * H[0] - partial_f[1] * H[1]) / detH;
        incr_y = (-2 * (x1 - x0) * H[0] + partial_f[0] * H[1]) / detH;

        /* Finally compute x_{k+1} and swap */
        x2 = x1 - incr_x;
        y2 = y1 - incr_y;

        x1 = x2;
        y1 = y2;

        n++;
    }

    if (n == MAX_ITERATIONS) exit_code = 2;

out:
    free(H);
    free(partial_f);

    result[0] = x1;
    result[1] = y1;

    return 0;
}

int main() {
    FILE *output;

    double h = H_STEP;
    int i, exit_code;

    /* Temp storage for the current point */
    double *point, *tgt;

    /* Store last iteration derivative to check we're going */
    /* in the right direction */
    double *tgt0;

    /* Time benchmarking */
    clock_t start_time, end_time;
    double elapsed_time;

    point = calloc(2, sizeof(double));
    if (point == NULL) return 2;

    tgt = calloc(2, sizeof(double));
    if (tgt == NULL) return 2;

    tgt0 = calloc(2, sizeof(double));
    if (tgt0 == NULL) return 2;

    output = fopen("results.txt", "w");
    if (output == NULL) return 3;

    start_time = clock();
    newton_zero_x(&point[0]);

    /* Save the first point */
    fprintf(output, "%.12lf %.12lf\n", point[0], point[1]);

    for (i = 0; i < PLOT_STEPS; i++) {
        exit_code = tangent(point, tgt);
        if (exit_code != 0) {
            fprintf(stderr, "Iter #%d: Tangent got division by 0 at (%f, %f)\n",
                    i, point[0], point[1]);
            break;
        }

        if (i > 0 && dot(tgt, tgt0) < 0) {
            tgt[0] *= -1;
            tgt[1] *= -1;
        }

        tgt0[0] = tgt[0];
        tgt0[1] = tgt[1];

        /* Newton method iteration to find the next point */
        exit_code = newton_iterator(point[0], point[1], point[0] + h * tgt[0],
                                    point[1] + h * tgt[1], h, point);

        if (exit_code != 0) {
            fprintf(stderr, "Iter #%d: Newton iterator returned code %d\n", i,
                    exit_code);
            break;
        }

        if (isnan(point[0]) || isnan(point[1])) {
            fprintf(stderr, "Iter #%d: Got NaN in current point: (%f, %f)\n", i,
                    point[0], point[1]);
            break;
        }

        if (fabs(f(point[0], point[1])) >= TOLERANCE) {
            fprintf(stderr, "Iter #%d: Tolerance exceeded (%2.16f)\n", i,
                    fabs(f(point[0], point[1])));
            break;
        }

        /* Save the point to a file */
        fprintf(output, "%.12lf %.12lf\n", point[0], point[1]);
    }

    end_time = clock();
    elapsed_time = ((double)end_time - (double)start_time) / CLOCKS_PER_SEC;
    fprintf(stderr, "Elapsed time: %2.4f s.\n", elapsed_time);

    /* Close output file and free memory */
    fclose(output);

    free(point);
    free(tgt);
    free(tgt0);

    return 0;
}
