
#Date: 18.3.24
#Group members: Andrey Romanchuk ID;Shahar Dahan ID;Maya Rozenberg ID
#Source Git: https://github.com/lihiSabag/Numerical-Analysis-2023.git
#Private Git:https://github.com/ShaharD25/Numerical-analysis-test-2.git
#Name: shahar dahan

import numpy as np

"""
Receives 3 parameters:
    1.  a - start value.
    2.  b - end  value. 
    3.  err - value of tolerable error

Returns variables:
    1.  S - The minimum number of iterations required to reach the desired accuracy
"""


def max_steps(a, b, err):
    s = int(np.floor(- np.log2(err / (b - a)) / np.log2(2) - 1))
    return s



def bisection_method(f, a, b, tol=1e-6):
    # if np.sign(a) == np.sign(b):
    #     raise Exception("The scalars a and b do not bound a root")
    c, k = 0, 0
    steps = max_steps(a, b, tol)  # calculate the max steps possible

    #print("{:<10} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}".format("Iteration", "a", "b", "f(a)", "f(b)", "c", "f(c)"))

    # while the diff af a&b is not smaller than tol, and k is not greater than the max possible steps
    while abs(b - a) > tol and k <= steps:
        c = (a + b) / 2  # Calculation of the middle value

        if f(c) == 0:
            return c  # Procedure completed successfully

        if f(c) * f(a) < 0:  # if sign changed between steps
            b = c  # move forward
        else:
            a = c  # move backward

        #print("{:<10} {:<15.6f} {:<15.6f} {:<15.6f} {:<15.6f} {:<15.6f} {:<15.6f}".format(k, a, b, f(a), f(b), c, f(c)))
        k += 1

    return c  # return the current root


def find_all_roots(f, a, b, tol=1e-6):
    roots = []
    x = np.linspace(a, b, 1000)  # Divide the interval into smaller sub-intervals

    for i in range(len(x) - 1):
        if np.sign(f(x[i])) != np.sign(f(x[i + 1])):
            root = bisection_method(f, x[i], x[i + 1], tol)
            roots.append(root)

    return roots


if __name__ == '__main__':
    f = lambda x: (4*x**2 - 8*x**2 ) / (x**2 + 4)

    # Adjust the interval to avoid the singularity
    a = -2
    b = 2

    roots = find_all_roots(f, a, b)
    print( f"\nThe equation f(x) has approximate roots at {roots}")