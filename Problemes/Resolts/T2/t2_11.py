import numpy as np


def compute_gershgorin(A):
    vap, err = np.inf, None

    for i, row in enumerate(A):
        v, r = row[i], sum(row) - row[i]

        # Compute the Gershgorin disks (intervals as in R)
        interval = np.array([r, -r]) if r < 0 else np.array([-r, r])
        interval += v

        if v < vap:
            vap = v
            err = interval

        print(f"Lambda({i}) in {interval}\t{v:2.4f} +- {np.abs(r):2.4f}")

    return vap, err


A = np.array([
    [2.14, -0.10, 0.00],
    [-0.10, 4.34, 0.20],
    [0.00, 0.20, 4.48],
])

print("(b)")
# Smallest is the first
compute_gershgorin(A)

print("\n(c)")
k, alpha = 0, 0.05

# Compute A_{k,alpha}
A[k, :] *= alpha
A[:, k] /= alpha

vap, err = compute_gershgorin(A)
best = min(np.linalg.eigvals(A))

print(f"\nReal min VAP = {best:2.6f}")

is_inside = 'Yes' if best >= err[0] and best <= err[1] else 'No'
print(f"Is in the interval? {is_inside}")

print(f"The error is {np.abs(best - vap):2.6f}")

print("\n(d)")
k, alpha = 1, 0.1

# Compute A_{k,alpha}
A[k, :] *= alpha
A[:, k] /= alpha

vap, err = compute_gershgorin(A)
vaps = np.linalg.eigvals(A)

print(f"\nReal VAPs = {vaps}")
