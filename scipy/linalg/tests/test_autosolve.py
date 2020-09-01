# coverage test for scipy autosolve ... about 99%
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import hilbert
from autosol import (autosolve, autosolve_nonneg,
                     autosolve_rising,two_norm)
# when autosol is online...
# from scipy.autosol import (autosolve, autosolve_nonneg,
#                   autosolve_rising, two_norm)


# inconsistent arrays .. these crash the run
# A=np.zeros((3,3))
# b=np.zeros(5)
# x=autosolve(A,b)[0]
# x=autosolve_nonneg(A,b)[0]
# x=autosolve_rising(A,b)[0]

# 0-D
print()
print("Test 1")
print("call autosolve with 0-dimensions")
x = autosolve(2.0, 2.0)[0]
print(x)
print("call autosolve_nonneg with 0-dimensions")
x = autosolve_nonneg(2.0, -2.0)[0]
print(x)
print("call autosolve_rising with 0-dimensions")
x = autosolve_rising(2.0, 2.0)[0]
print(x)

# 1 row
print()
print("Test 2")
print("test all three solvers for zero-dimensional matrix")
print("call autosolve with 1 row")
x = autosolve([[1.0, 2.0, 3.0]], 12.0)[0]
print(x)
print("call autosolve_nonneg with 1 row")
x = autosolve_nonneg([[1.0, 2.0, 3.0]], -12.0)[0]
print(x)
print("call autosolve_rising with 1 row")
x = autosolve_rising(2.0, -2.0)[0]
print(x)

# 1 column
print()
print("Test 3")
print("call autosolve with 1-column")
x = autosolve([[1.0], [2.0], [3.0]], [[3.0], [6.0], [9.0]])[0]
print(x)
print("call autosolveNN with 1-column")
x = autosolve_nonneg([[1.0], [2.0], [3.0]], [[3.0], [6.0], [9.0]])[0]
x = autosolve_nonneg([[1.0], [2.0], [3.0]], [[3.0], [6.0], [-9.0]])[0]
print(x)
print("call autosolve_rising with 1 column")
x = autosolve_rising([[1.0], [2.0], [3.0]], [[3.0], [6.0], [9.0]])[0]
print(x)

print()
print("Test 4")
print("six tests to check all zeros input")
A = np.zeros((3, 3))
A1 = np.zeros((3, 3))
A1[0, 1] = 1.0
b = np.zeros(3)
b1 = np.zeros(3)
b1[1] = 1.0
x = autosolve(A, b1)
print(x)
x = autosolve_nonneg(A, b1)
print(x)
x = autosolve_rising(A, b1)
print(x)
x = autosolve(A1, b)
print(x)
x = autosolve_nonneg(A1, b)
print(x)
x = autosolve_rising(A1, b)
print(x)

print()
print("Test 5")
print("test with real zero singular value.")
A = np.zeros((3, 3))
A[0, 0] = 1
b = [1.0, 0.0, 0.0]
x = autosolve(A, b)[0]
print(x)

print()
print("Test 6")
print("see if nonneg handles a problem that is trivially nonneg.")
A = np.ones((2, 2))
A[1, 1] = 2
x = np.ones(2)
x[1] = 2
b = A @ x
x = autosolve_nonneg(A, b)[0]
print(x)
x = autosolve_rising(A, b)[0]
print(x)

print()
print("Test 7")
print("demo autosolve on three different matix sizes")
for m in (5, 7, 25):
    n = m
    # plot axis
    ax = np.zeros(n)
    for i in range(0, n):
        ax[i] = float(i) / (float(n) - 1.0)
    A = hilbert(n)
    truex = np.zeros(n)
    if n < 4:
        for i in range(0, n):
            truex[i] = float(i + 1)
    else:
        cut = int(float(n) * 1.26)
        for i in range(0, n):
            truex[i] = float(min(i + 1, cut - i))
    b = np.matmul(A, truex)
    # add extra noise
    for i in range(0, m):
        b[i] += 0.00005 * np.sin(1.2 + 4.0 * i)
    x, ur, sigma, lambdah = autosolve(A, b)
    print("autosolve solution for n= ", n)
    print(x)
    #plt.plot(truex)
    #plt.show()
    #plt.plot(x)
    #plt.show()


print()
print("Test 8")
print("demo nonnegative with inverse heat problem.")
for m in (5, 15):
    n = m
    A = np.zeros((m, n))
    b = np.zeros(m)
    truex = np.zeros(n)
    for i in range(0, m):
        s = 1.5 * float(i) / float(m - 1)
        for j in range(0, n):
            t = 1.5 * float(j) / float(n - 1)
            A[i, j] = np.exp(-(s - t) * (s - t))
    for i in range(0, n):
        truex[i] = min(i, n + 1 - i) - 1
    b = A @ truex
    for i in range(0, m):
        b[i] += 0.00005 * np.sin(1.4 + 4.2 * i)
    x = autosolve_nonneg(A, b)[0]
    print("autosolve_nonneg solution for n= ", n)
    print(x)
    #plt.plot(truex)
    #plt.show()
    #plt.plot(x)
    #plt.show()


print()
print("Test 9")
print("demo autosolve_rising with wiggly solution")

for m in (15, 25):
    n = m
    A = np.zeros((m, n))
    b = np.zeros(m)
    truex = np.zeros(n)
    for i in range(0, m):
        s = 1.5 * float(i) / float(m - 1)
        for j in range(0, n):
            t = 1.5 * float(j) / float(n - 1)
            A[i, j] = np.exp(-(s - t) * (s - t))
    for i in range(0, n):
        truex[i] = 0.75 * float(i) + np.sin(float(i))
    b = A @ truex
    x = autosolve_rising(A, b)[0]
    print("autosolve_rising solution")
    print(x)
    #plt.plot(truex)
    #plt.show()
    #plt.plot(x)
    #plt.show()

print()
print("Test 10")
print("test two-norm for tiny arrays")
b=np.ones(2)
print(two_norm(b))
b=np.ones(1)
print(two_norm(b))
b=np.ones(0)
print(two_norm(b))
#2345678901234567890123456789012345678901234567890123456789012345678901234567890
'''
Test 1
call autosolve with 0-dimensions
[1.]
call autosolve_nonneg with 0-dimensions
[0.]
call autosolve_rising with 0-dimensions
1.0

Test 2
test all three solvers for zero-dimensional matrix
call autosolve with 1 row
[0.85714286 1.71428571 2.57142857]
call autosolve_nonneg with 1 row
[0. 0. 0.]
call autosolve_rising with 1 row
-1.0

Test 3
call autosolve with 1-column
[[3.]]
call autosolveNN with 1-column
[0.]
call autosolve_rising with 1 column
[3.]

Test 4
six tests to check all zeros input
(array([0., 0., 0.]), 0, 0.0, 0.0)
(array([0., 0., 0.]), 0, 0.0, 0.0)
(array([0., 0., 0.]), 0, 0.0, 0.0)
(array([0., 0., 0.]), 0, 0.0, 0.0)
(array([0., 0., 0.]), 0, 0.0, 0.0)
(array([0., 0., 0.]), 0, 0.0, 0.0)

Test 5
test with real zero singular value.
[1. 0. 0.]

Test 6
see if nonneg handles a problem that is trivially nonneg.
1.0000000000000009
1.0000000000000009

Test 7
demo autosolve on three different matix sizes
autosolve solution for n=  5
[1.0022663  1.96715913 3.1323387  2.80343516 2.09613139]
autosolve solution for n=  7
[1.01023617 1.8208179  3.56546253 3.79063197 3.44517575 2.94513422
 2.44350035]
autosolve solution for n=  25
[ 1.02259448  1.58214106  4.38265186  3.64309312  3.77144851  4.95939886
  6.67073887  8.46653312 10.09362302 11.43500288 12.45541031 13.16363538
 13.58962851 13.77149077 13.7484734  13.55746903 13.23146872 12.79908993
 12.28466072 11.70856972 11.08772108 10.43600845  9.7647646   9.08316711
  8.39859361]

Test 8
demo nonnegative with inverse heat problem.
autosolve_nonneg solution for n=  5
[0.         0.         0.         1.67299927 1.90606466]
autosolve_nonneg solution for n=  15
[0.         0.         0.         0.         2.09213789 7.18369773
 7.57157313 5.98075872 4.4569618  4.01516893 4.59550375 5.30308821
 4.83678217 1.97247424 0.        ]

Test 9
demo autosolve_rising with wiggly solution
autosolve_rising solution
[-1.98734984e-03  1.60382527e+00  2.38381726e+00  2.38381726e+00
  2.38381726e+00  2.47464790e+00  4.58972561e+00  5.66370084e+00
  7.08571870e+00  7.08571870e+00  7.08571870e+00  7.11709151e+00
  8.53944877e+00  1.01465884e+01  1.14937508e+01]
autosolve_rising solution
[7.01073301e-03 1.57132601e+00 2.40680970e+00 2.40680970e+00
 2.40680970e+00 2.40680970e+00 4.36897617e+00 6.34754874e+00
 6.34754874e+00 7.32014107e+00 7.32014107e+00 7.32014107e+00
 7.32014107e+00 1.13513111e+01 1.13513111e+01 1.16445585e+01
 1.16445585e+01 1.16445585e+01 1.35049269e+01 1.35049269e+01
 1.64484244e+01 1.64484244e+01 1.64484244e+01 1.64484244e+01
 1.70841648e+01]

Test 10
test two-norm for tiny arrays
1.4142135623730951
1.0
0.0
'''
