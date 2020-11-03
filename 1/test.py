import numpy as np

N = 1000
ERROR = 1e-12
MAX_ITER = 1000


def jacobi(A, b, x_init, epsilon=ERROR, max_iterations=MAX_ITER):
    D = np.diag(np.diag(A))
    LU = A - D
    x = x_init
    for i in range(max_iterations):
        D_inv = np.diag(1 / np.diag(D))
        x_new = np.dot(D_inv, b - np.dot(LU, x))
        if np.linalg.norm(x_new - x) < epsilon / 4.:
            return x_new
        x = x_new
    return x


def get_matrix_value(i, j, n=N):
    if i == j:
        return 5

    if j + 3 == i or j - 3 == i:
        return 2

    if (i == 0 and j == n - 1) or (i == n - 1 and j == 0):
        return 1

    return 0


def main():
    from time import time

    t = time()

    A = np.fromiter((
        get_matrix_value(i, j) for j in range(N)
        for i in range(N)), float).reshape((N, N))
    b = np.fromiter(((i + 1) / N for i in range(N)), float)

    t = time() - t
    print(f"Init:  {t:2.4f} s.")

    t = time()

    # x_init = np.zeros(N)
    # x = jacobi(A, b, x_init)
    x = np.linalg.solve(A, b)

    t = time() - t
    print(f"Solve: {t:2.4f} s.")

    # Assert Numpy error is actually ok
    # assert np.linalg.norm(np.dot(A, x) - b) < ERROR

    file_names = [
        "./Jacobi_YaremaTaras.txt",
        "./Gauss-Seidel_YaremaTaras.txt"
    ]

    for file_name in file_names:
        with open(file_name, "r") as f:
            target = np.array(
                [np.float(e) for e in f.readlines()[0].split(',')]
            )

            diff = np.linalg.norm(target - x)

            # Check for the error we got with C
            try:
                assert diff < ERROR
            except AssertionError:
                print(f"Large error ({file_name}) = {diff}")
                exit(1)

    print("OK")
    exit(0)


if __name__ == "__main__":
    main()
