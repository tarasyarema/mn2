import numpy as np
from numpy.linalg import norm as norm
from numpy.linalg import inv as inv


def newton(x: np.array, df, f, error=1e-8):
    its = 0

    while True:
        y = x - inv(df(x)) @ f(x)
        its += 1

        if norm(y - x, 2) < error or norm(f(y), 2) < error:
            break

        x = y

    print("NEWTON")
    print(f"Took {its} iterations")
    print(f"Sol: {x}")


def quasi_newton(x: np.array, df, f, error=1e-8):
    its = 0
    _df = inv(df(x))

    while True:
        y = x - _df @ f(x)
        its += 1

        if norm(y - x, 2) < error or norm(f(y), 2) < error:
            break

        x = y

    print("QUASI NEWTON")
    print(f"Took {its} iterations")
    print(f"Sol: {x}")


def f(arr: np.array):
    x, y = arr
    return np.array([
        x**3 + y**2 + x + y - 2,
        x*y**2 + x**2 + x - y
    ])


def df(arr: np.array):
    x, y = arr
    return np.array([
        [3*x**2 + 1, 2*y + 1],
        [y**2 + 2*x + 1, 2*x*y - 1],
    ])


x = np.array([0.4, 0.8])

newton(x, df, f)
quasi_newton(x, df, f)
