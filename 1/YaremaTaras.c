#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Handles if the execution mode is
 * for submitting */
#define SUBMIT 0

#if SUBMIT
#define N 1000000
#else
/* This is used for debugging and testing purposes */
#define N 1000
#endif

/* The tolerance, or acceptable error */
#define ERROR (double)1e-12

/* Max iterations to mark proccess as not convergent */
#define ITER_MAX 1000

/* Output files as constants */
const char *jacobi_file = "Jacobi_YaremaTaras.txt";
const char *gauss_seidel_file = "Gauss-Seidel_YaremaTaras.txt";
const char *sor_file = "SOR_YaremaTaras.txt";

/* Prints the resume of a function execution */
void print_resume(const char *n, int its, double err, double t) {
#if SUBMIT
    printf("Algoritme: %s Iterats: %d Error: %.8e\n", n, its, err);
#else
    printf("%15s | its = %4d | err = %.8e | time = %2.6f s\n", n, its, err, t);
#endif
}

/* Prints the solution vector in the correct format */
int print_solution(const char *name, double *sol) {
    FILE *f = (FILE *)fopen(name, "w");
    int i;

    /* Try to open the file */
    if (f == NULL) {
        printf("could not open file '%s'\n", name);
        return 1;
    }

    /* Print in the correct format */
#if SUBMIT
    for (i = 0; i < N; i += 1000) {
#else
    for (i = 0; i < N; ++i) {
#endif
        fprintf(f, "%s%2.14f", i == 0 ? "" : ",", sol[i]);
    }

    fprintf(f, "\n");

    /* Close file */
    fclose(f);

    return 0;
}

/* Compute the element A[i][j] */
double get_matrix_value(int i, int j) {
    /* We assue that you give values in the correct ranges,
     * i.e. 0 <= i < N and 0 <= j < N. */
    if (i == j) return 5;
    if (j + 3 == i || j - 3 == i) return 2;
    if ((i == 0 && j == N - 1) || (i == N - 1 && j == 0)) return 1;
    return 0;
}

/* Compute the element b[i] */
double get_independent_value(int i) {
    /* We assue that you give value in the correct range,
     * i.e. 0 <= i < N. */
    return (double)(i + 1) / (double)N;
}

/* Compute the element i-th element of D^-1 * b, c[i] */
double get_c_value(int i) {
    /* We assue that you give value in the correct range,
     * i.e. 0 <= i < N. */
    return get_independent_value(i) / get_matrix_value(i, i);
}

/* Compute the element B[i][j] of the iteration matrix */
double get_bj_matrix_value(int i, int j) {
    /* We assue that you give values in the correct ranges,
     * i.e. 0 <= i < N and 0 <= j < N. */
    if (j + 3 == i || j - 3 == i) return 2. / 5.;
    if ((i == 0 && j == N - 1) || (i == N - 1 && j == 0)) return 1. / 5.;
    return 0;
}

/* Computes the norm of the matrix given by the function getter.
 * Implemented getters:
 *      get_matrix_value (original matrix A),
 *      get_bj_matrix_value (iterarion matrix B)
 */
double matrix_norm_inf(double (*getter)(int, int)) {
    double best = 0, partial;
    int i;

    for (i = 0; i < N; ++i) {
        /* Compute the partial sum via columns.
         * Allways get the diagonal value and the first upper diag. value */
        partial = (*getter)(i, i) + (*getter)(0, 3);

        /* Get the left-bot and right-top values */
        if (i == 0) partial += (*getter)(i, N - 1);
        if (i == N - 1) partial += (*getter)(0, i);

        /* Get the additional lower/upper diagonal value if needed */
        if (i >= 3 && i + 3 < N) partial += (*getter)(0, 3);

        /* Get absolute value */
        partial = fabs(partial);

        /* Update best */
        best = partial > best ? partial : best;
    }

    return best;
}

/* Computes the norm of the vector v */
double vector_norm(double *v) {
    double norm = 0;
    int i = 0;
    for (i = 0; i < N; ++i) norm += (v[i] * v[i]);
    return sqrt(norm);
}

/* Computes the ||δ(k)||_2 = ||actual - last|| = ||x(k+1) - x(k)|| */
double delta_norm(double *actual, double *last) {
    double norm = 0, tmp = 0;
    int i = 0;
    for (i = 0; i < N; ++i) {
        tmp = actual[i] - last[i];
        norm += (tmp * tmp);
    }
    return sqrt(norm);
}

/* Computes the ||δ(k)||_inf = ||actual - last|| = ||x(k+1) - x(k)|| */
double delta_norm_inf(double *actual, double *last) {
    double norm = 0, tmp;
    int i = 0;

    for (i = 0; i < N; ++i) {
        tmp = fabs(actual[i] - last[i]);
        if (tmp > norm) norm = tmp;
    }

    return norm;
}

/* Utility function two swap to pointers of type *double */
void swap(double **a, double **b) {
    double *temp = *a;
    *a = *b;
    *b = temp;
}

/* Computes the Jacobi method for the given stationary system.
 * Total allocs: 2 * N * 8 bytes, i.e. for N = 10^6 we need 0.16 megabytes. */
int jacobi() {
    double delta, elapsed;
    int i, iterations = 0;

    /* Benchmarking purposes */
    clock_t t;

    /* Init the solution vectors */
    double **sol = malloc(2 * sizeof(double *));
    sol[0] = (double *)calloc(N, sizeof(double));
    sol[1] = (double *)calloc(N, sizeof(double));

    t = clock();

    while (iterations < ITER_MAX) {
        /* This method is optiomal because we do not store any matrix.
         * So to compute every element of the solution vector we do so in O(N)
         * time, no O(n^2) as we would with a matrix. */
        for (i = 0; i < N; ++i) {
            sol[1][i] = get_c_value(i);

            /* First row needs to sum the right-top corner */
            if (i == 0)
                sol[1][i] -= get_bj_matrix_value(i, N - 1) * sol[0][N - 1];

            /* Last row needs to sum the left-bot corner */
            if (i == N - 1) sol[1][i] -= get_bj_matrix_value(0, i) * sol[0][0];

            /* Lower diagonal element */
            if (i >= 3)
                sol[1][i] -= get_bj_matrix_value(i, i - 3) * sol[0][i - 3];

            /* Upper diagonal element */
            if (N - i > 3)
                sol[1][i] -= get_bj_matrix_value(0, 3) * sol[0][i + 3];
        }

        /* Compute the actual absolute error */
        delta = delta_norm_inf(sol[1], sol[0]);

        /* Swap the solutions */
        swap(&sol[0], &sol[1]);
        /* for (i = 0; i < N; ++i) sol[1][i] = sol[0][i]; */

        /* Check for stop criterion, note the 4 dividing.
         * Its computed via the absolute error stop criterion with the
         * norm of the matrix B = 0.80 */
        if (iterations > 2 && delta <= ERROR / 4.) {
            break;
        }

        ++iterations;
    }

    /* Stdout debugging */
    elapsed = (double)(clock() - t) / (double)CLOCKS_PER_SEC;
    print_resume("jacobi", iterations, delta, elapsed);

    /* Print solution to correct file */
    print_solution(jacobi_file, sol[0]);

    /* Free everything */
    free(sol[0]);
    free(sol[1]);
    free(sol);

    return iterations;
}

/* Computes the Gauss-Seidel method for the given stationary system.
 * Total allocs: 2 * N * 8 bytes, i.e. for N = 10^6 we need 0.16 megabytes. */
int gauss_seidel() {
    double delta, elapsed;
    int i, iterations = 0;

    /* Benchmarking purposes */
    clock_t t;

    /* Init the solution vectors */
    double **sol = malloc(2 * sizeof(double *));

    /* The k-th x */
    sol[0] = (double *)calloc(N, sizeof(double));

    /* The (k+1)-th x */
    sol[1] = (double *)calloc(N, sizeof(double));

    t = clock();

    while (iterations < ITER_MAX) {
        /* As in Jacobi, this method is also O(n) */
        for (i = 0; i < N; ++i) {
            sol[1][i] = get_c_value(i);

            /* First row needs to sum the right-top corner */
            if (i == 0)
                sol[1][i] -= get_bj_matrix_value(i, N - 1) * sol[0][N - 1];

            /* Last row needs to sum the left-bot corner */
            if (i == N - 1) sol[1][i] -= get_bj_matrix_value(0, i) * sol[1][0];

            /* Lower diagonal element */
            if (i >= 3)
                sol[1][i] -= get_bj_matrix_value(i, i - 3) * sol[1][i - 3];

            /* Upper diagonal element */
            if (N - i > 3)
                sol[1][i] -= get_bj_matrix_value(i, i + 3) * sol[0][i + 3];
        }

        /* Compute the actual absolute error */
        delta = delta_norm_inf(sol[1], sol[0]);

        /* Swap the solutions */
        swap(&sol[0], &sol[1]);
        /* for (i = 0; i < N; ++i) sol[1][i] = sol[0][i]; */

        /* Check for stop criterion, note the 2 dividing.
         * Its computed via the absolute error stop criterion with the
         * norm of the matrix B <= 2/3 */
        if (iterations > 2 && delta <= ERROR / 4.) {
            break;
        }

        ++iterations;
    }

    /* Stdout debugging */
    elapsed = (double)(clock() - t) / (double)CLOCKS_PER_SEC;
    print_resume("gauss-seidel", iterations, delta, elapsed);

    /* Print solution to correct file */
    print_solution(gauss_seidel_file, sol[0]);

    /* Free everything */
    free(sol[0]);
    free(sol[1]);
    free(sol);

    return iterations;
}

int sor(double w, int least_iters, int want_print) {
    double delta, elapsed;
    int i, iterations = 0;
    char name_buffer[16];

    /* Benchmarking purposes */
    clock_t t;

    /* Init the solution vectors */
    double **sol = malloc(2 * sizeof(double *));

    /* The k-th x */
    sol[0] = (double *)calloc(N, sizeof(double));

    /* The (k+1)-th x */
    sol[1] = (double *)calloc(N, sizeof(double));

    t = clock();

    while (iterations < ITER_MAX) {
        /* As in Jacobi, this method is also O(n) */
        for (i = 0; i < N; ++i) {
            sol[1][i] = get_c_value(i);

            /* First row needs to sum the right-top corner */
            if (i == 0)
                sol[1][i] -= get_bj_matrix_value(i, N - 1) * sol[0][N - 1];

            /* Last row needs to sum the left-bot corner */
            if (i == N - 1) sol[1][i] -= get_bj_matrix_value(0, i) * sol[1][0];

            /* Lower diagonal element */
            if (i >= 3)
                sol[1][i] -= get_bj_matrix_value(i, i - 3) * sol[1][i - 3];

            /* Upper diagonal element */
            if (N - i > 3)
                sol[1][i] -= get_bj_matrix_value(i, i + 3) * sol[0][i + 3];

            /* Compute right side */
            sol[1][i] = w * (sol[1][i] - sol[0][i]);

            /* Add left side */
            sol[1][i] += sol[0][i];
        }

        /* Compute the actual absolute error */
        delta = delta_norm_inf(sol[1], sol[0]);

        /* Swap the solutions */
        swap(&sol[0], &sol[1]);
        /* for (i = 0; i < N; ++i) sol[1][i] = sol[0][i]; */

        /* Check for stop criterion, note the 2 dividing.
         * Its computed via the absolute error stop criterion with the
         * norm of the matrix B <= 2/3 */
        if (iterations > 2 && delta <= ERROR / 4.) {
            break;
        }

        ++iterations;
    }

    /* Stdout debugging */
    elapsed = (double)(clock() - t) / (double)CLOCKS_PER_SEC;

    if (want_print) {
        /* Sprintf the sor name */
        sprintf(name_buffer, "sor(%1.6f)", w);
        print_resume(name_buffer, iterations, delta, elapsed);

        /* Only print solution when asked for */
        print_solution(sor_file, sol[0]);
    }

    /* Free everything */
    free(sol[0]);
    free(sol[1]);
    free(sol);

    return iterations;
}

int main() {
    /* Best computed SOR with step = 0.0001 is w = 1.220600,
     * elapsing around 100 seconds */
    double w, w_best, step = 0.005, elapsed;
    int sor_actual, sor_best = ITER_MAX, iters = 0;
    clock_t t;

    jacobi();
    gauss_seidel();

    t = clock();

    for (w = 0.1; w < 2; w += step) {
        sor_actual = sor(w, sor_best, 0);

        if (sor_actual < sor_best) {
            sor_best = sor_actual;
            w_best = w;
        }

        ++iters;
    }

    elapsed = (double)(clock() - t) / (double)CLOCKS_PER_SEC;

    sor(w_best, ITER_MAX, 1);
    printf("%15s | %d steps computed in %2.4f s\n", "sor", iters, elapsed);

    return 0;
}
