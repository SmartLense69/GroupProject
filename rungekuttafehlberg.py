import numpy as np
import matplotlib.pyplot as plt
import diffequation as de


def diff1(inXi, inPhi, inTheta, n=1):
    if inXi == 0:
        return 0
    return - inTheta**n - (2/inXi) * inPhi


def diff2(inPhi):
    return inPhi


A0 = np.asarray([
    [None, None, None, None, None, None, None],
    [None, None, None, None, None, None, None],
    [None, (1 / 5), None, None, None, None, None],
    [None, (3 / 40), (9 / 40), None, None, None, None],
    [None, (44 / 45), (-56 / 15), (32 / 9), None, None, None],
    [None, (19372 / 6561), (-25360 / 2187), (64448 / 6561), (-212 / 729), None, None],
    [None, (9017 / 3168), (-355 / 33), (46732 / 5247), (49 / 176), (-5103 / 18656), None],
    [None, (35 / 384), 0, (500 / 1113), (125 / 192), (-2187 / 6784), (11 / 84)]
], dtype=np.float64)

A = np.asarray([
    [None, None, None, None, None, None, None],
    [None, None, None, None, None, None, None],
    [None, (1 / 4), None, None, None, None, None],
    [None, (3 / 32), (9 / 32), None, None, None, None],
    [None, (1932 / 2197), (-7200 / 2197), (7296 / 2197), None, None, None],
    [None, (439 / 216), (-8), (3680 / 513), (-845 / 4104), None, None],
    [None, (-8 / 27), (2), (-3544 / 2565), (1859 / 4104), (-11 / 40), None],
], dtype=np.float64)
"""
"a" Coefficients for Runge-Kutta-Fehlberg.
Index the array like the index mentioned in the math scripts, i.e.
first element starts at 1.
This was made to make it easier to copy the algorithm.
"""

B0 = np.asarray(
    [None, (35 / 384), 0, (500 / 1113), (125 / 192), (-2187 / 6784), (11 / 84), 0], dtype=np.float64)

B = np.asarray(
    [None, (16 / 135), 0, (6656 / 12825), (28561 / 56430), (-9 / 50), (2 / 55)], dtype=np.float64)
"""
"b" Coefficients for Runge-Kutta-Fehlberg.
"""

Bm = np.asarray(
    [None, (5179 / 57600), 0, (7517 / 16695), (393 / 640), (-92097 / 339200), (187 / 2100), (1 / 40)], dtype=np.float64)
"""
"b*" Coefficients for Runge-Kutta-Fehlberg.
"""

C0 = np.asarray(
    [None, 0, (1 / 5), (3 / 10), (4 / 5), (8 / 9), 1, 1], dtype=np.float64)

C = np.asarray(
    [None, 0, (1 / 4), (3 / 8), (12 / 13), 1, (1 / 2)], dtype=np.float64)
"""
"c" Coefficients for Runge-Kutta-Fehlberg.
"""


def rkf(diffEq1: de.DifferentialEquation, diffEq2: de.DifferentialEquation, stepSize, maxTime=1, coeff=1):
    """
    Solves a differential equation system.
    This function only applies for Lane-Emden currently.
    :param diffEq1: First decoupled diff. equation (-theta^4+2phi/xi)
    :param diffEq2: Second decoupled diff. equation (phi)
    :param stepSize: Size of the step, often mentioned in the scripts as h
    :param maxTime: The max time the differential equation should be solved to.
    :param coeff: Additional parameter for Lane-Emden
    :return: xi, phi and theta of the set of differential equations
    """

    xi = np.arange(0, maxTime + stepSize, stepSize)
    n = np.shape(xi)[0]
    phi = np.zeros(n)
    theta = np.zeros(n)
    phi[0] = diffEq1.iniC
    theta[0] = diffEq2.iniC

    for i in range(1, n):

        # Solve the first differential equation
        # Keep theta constant

        j1 = stepSize * diffEq1.func(xi[i - 1],
                                     phi[i - 1],
                                     theta[i - 1], coeff)
        k1 = stepSize * (phi[i - 1])

        j2 = stepSize * diffEq1.func(xi[i - 1] + C[2] * stepSize,
                                     phi[i - 1] + A[2, 1] * j1,
                                     theta[i - 1] + A[2, 1] * k1, coeff)
        k2 = stepSize * (phi[i - 1] + A[2, 1] * j1)

        j3 = stepSize * diffEq1.func(xi[i - 1] + C[3] * stepSize,
                                     phi[i - 1] + A[3, 1] * j1 + A[3, 2] * j2,
                                     theta[i - 1] + A[3, 1] * k1 + A[3, 2] * k2, coeff)
        k3 = stepSize * (phi[i - 1] + A[3, 1] * j1 + A[3, 2] * j2)

        j4 = stepSize * diffEq1.func(xi[i - 1] + C[4] * stepSize,
                                     phi[i - 1] + A[4, 1] * j1 + A[4, 2] * j2 + A[4, 3] * j3,
                                     theta[i - 1] + A[4, 1] * k1 + A[4, 2] * k2 + A[4, 3] * k3, coeff)
        k4 = stepSize * (phi[i - 1] + A[4, 1] * j1 + A[4, 2] * j2 + A[4, 3] * j3)

        j5 = stepSize * diffEq1.func(xi[i - 1] + C[5] * stepSize,
                                     phi[i - 1] + A[5, 1] * j1 + A[5, 2] * j2 + A[5, 3] * j3 + A[5, 4] * j4,
                                     theta[i - 1] + A[5, 1] * k1 + A[5, 2] * k2 + A[5, 3] * k3 + A[5, 4] * k4, coeff)
        k5 = stepSize * (phi[i - 1] + A[5, 1] * j1 + A[5, 2] * j2 + A[5, 3] * j3 + A[5, 4] * j4)

        j6 = stepSize * diffEq1.func(xi[i - 1] + C[6] * stepSize,
                                     phi[i - 1] + A[6, 1] * j1 + A[6, 2] * j2 + A[6, 3] * j3 + A[6, 4] * j4 + A[6, 5] * j5,
                                     theta[i - 1] + A[6, 1] * k1 + A[6, 2] * k2 + A[6, 3] * k3 + A[6, 4] * k4 + A[6, 5] * k5, coeff)
        k6 = stepSize * (phi[i - 1] + A[6, 1] * j1 + A[6, 2] * j2 + A[6, 3] * j3 + A[6, 4] * j4 + A[6, 5] * j5)

        phi[i] = phi[i - 1] + B[1] * j1 + B[2] * j2 + B[3] * j3 + B[4] * j4 + B[5] * j5 + B[6] * j6
        theta[i] = theta[i - 1] + B[1] * k1 + B[2] * k2 + B[3] * k3 + B[4] * k4 + B[5] * k5 + B[6] * k6

    return xi, phi, theta


def plot(n):

    # Define two differential equations, which will
    # represent the two decoupled Lane-Emden eqs.

    diffEq1 = de.DifferentialEquation(diff1, 0)
    diffEq2 = de.DifferentialEquation(diff2, 1)

    # Solve with RKF
    xi, phi, theta = rkf(diffEq1, diffEq2, 0.001, 3, n)

    # Plot
    nLabel = "$n=" + str(n) + "$ with RKF"
    plt.plot(xi, theta, label=nLabel)