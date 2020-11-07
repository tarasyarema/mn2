import numpy as np


# From Wikipedia
def power_iteration(A, num_simulations: int):
    b_k = np.random.rand(A.shape[1])

    for _ in range(num_simulations):
        b_k1 = np.dot(A, b_k)
        b_k1_norm = np.linalg.norm(b_k1)
        b_k = b_k1 / b_k1_norm

    return b_k


A = np.array([
    [+0, -1, -1],
    [+1, +1, -1],
    [+1, +1, +1],
])

vep = power_iteration(A, 100)

# Use the Rayleigh quotient
vap = np.dot(np.transpose(vep), np.dot(A, vep))
vap /= np.dot(np.transpose(vep), vep)

print("(a)")
print(f"vap = {vap}")
print(f"vep = {vep}")
print("not convergent...")

vep = power_iteration(np.linalg.inv(A), 100)

# Use the Rayleigh quotient
vap = np.dot(np.transpose(vep), np.dot(A, vep))
vap /= np.dot(np.transpose(vep), vep)

print("(b)")
print(f"vap_1 = {vap}")
print(f"vep_1 = {vep}")

det = np.linalg.det(A)
tr = sum(np.diag(A))
det /= vap
tr -= vap

print(f"vap_2 * vap_3 = {det}")
print(f"vap_2 + vap_3 = {tr}")

y = np.roots([1, -tr, det])
print(f"vap_2 = {y[0]}")
print(f"vap_3 = {tr - y[0]}")
