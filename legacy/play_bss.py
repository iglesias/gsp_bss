#!/usr/bin/python3

import cvxpy as cvx
import numpy as np
import tensorly

def PsiiTimeshi(Psi, h):
  return np.diag(np.matmul(Psi, h))

# Problem data.

# Shift operator.
S = np.array([[0, 1], [1, 0]])
[Delta, V] = np.linalg.eig(S)
U = np.linalg.inv(V)

num_filters = 2
# Number of coefficients in each filter.
L = np.array([2, 2])
assert(L.size == num_filters)
# Filter coefficients.
h = list()
for i in range(num_filters):
  h.append(np.random.randn(L[i]))

# Filters.
#TODO generalize.
H1 = h[0][0] + h[0][1]*S
H2 = h[1][0] + h[1][1]*S

Psi1 = np.vander(Lambda, len(h[0]), increasing=True)
Psi2 = np.vander(Lambda, len(h[1]), increasing=True)

# Input signals.
x1_truth = np.random.rand(2)
x2_truth = np.random.rand(2)

# Output signal.
y = np.dot(H1, x1_truth) + np.dot(H2, x2_truth)

PsiTimesh = np.concatenate((PsiiTimeshi(Psi1, h[0]), PsiiTimeshi(Psi2, h[1])), axis=1)
UTimesx = np.concatenate((np.dot(U, x1_truth), np.dot(U, x2_truth)))
ytilde = np.dot(PsiTimesh, UTimesx)

# Construct the problem.

x1 = cvx.Variable(2)
x2 = cvx.Variable(2)

# Using cvx.kron and reshape is only a fancy way to do the outer product.
Z1 = cvx.reshape(cvx.kron(h[0], x1), 2, 2)
Z2 = cvx.reshape(cvx.kron(h[1], x2), 2, 2)

Z = Z1+Z2

objective = cvx.Minimize(cvx.norm(Z, "nuc"))

A = tensorly.tenalg.khatri_rao((Psi1.T, U.T)).T

# Nonconvex warning.
# non_overlapping_constraint = np.inner(x1, x2.T) == 0
