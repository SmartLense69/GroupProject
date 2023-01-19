import numpy as np
import diffequation as DE


A = np.asarray([
    [None, None, None, None, None, None, None],
    [None, None, None, None, None, None, None],
    [None, (1 / 5), None, None, None, None, None],
    [None, (3 / 40), (9 / 40), None, None, None, None],
    [None, (44 / 45), (-56 / 15), (32 / 9), None, None, None],
    [None, (19372 / 6561), (-25360 / 2187), (64448 / 6561), (-212 / 729), None, None],
    [None, (9017 / 3168), (-355 / 33), (46732 / 5247), (49 / 176), (-5103 / 18656), None],
    [None, (35 / 384), 0, (500 / 1113), (125 / 192), (-2187 / 6784), (11 / 84)]
], dtype=np.float64)
"""
"a" Coefficients for Runge-Kutta-Fehlberg.
Index the array like the index mentioned in the math scripts, i.e.
first element starts at 1.
This was made to make it easier to copy the algorithm.
"""

B = np.asarray(
    [None, (35 / 384), 0, (500 / 1113), (125 / 192), (-2187 / 6784), (11 / 84), 0], dtype=np.float64)
"""
"b" Coefficients for Runge-Kutta-Fehlberg.
"""

Bm = np.asarray(
    [None, (5179 / 57600), 0, (7517 / 16695), (393 / 640), (-92097 / 339200), (187 / 2100), (1 / 40)], dtype=np.float64)
"""
"b*" Coefficients for Runge-Kutta-Fehlberg.
"""

C = np.asarray(
    [None, 0, (1 / 5), (3 / 10), (4 / 5), (8 / 9), 1, 1], dtype=np.float64)
"""
"c" Coefficients for Runge-Kutta-Fehlberg.
"""


def rkf(diffEq1: DE.DifferentialEquation, diffEq2: DE.DifferentialEquation, stepSize, maxTime=1, coeff=1):
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

        k1 = stepSize * diffEq1.func(xi[i], phi[i - 1], theta[i - 1], coeff)
        k2 = stepSize * diffEq1.func(xi[i] + C[2] * stepSize,
                                     phi[i - 1] + A[2, 1] * k1, theta[i - 1], coeff)
        k3 = stepSize * diffEq1.func(xi[i] + C[3] * stepSize,
                                     phi[i - 1] + A[3, 1] * k1 + A[3, 2] * k2, theta[i - 1], coeff)
        k4 = stepSize * diffEq1.func(xi[i] + C[4] * stepSize,
                                     phi[i - 1] + A[4, 1] * k1 + A[4, 2] * k2 + A[4, 3] * k3, theta[i - 1], coeff)
        k5 = stepSize * diffEq1.func(xi[i] + C[5] * stepSize,
                                     phi[i - 1] + A[5, 1] * k1 + A[5, 2] * k2 + A[5, 3] * k3 + A[5, 4] * k4,
                                     theta[i - 1], coeff)
        k6 = stepSize * diffEq1.func(xi[i] + C[6] * stepSize,
                                     phi[i - 1] + A[6, 1] * k1 + A[6, 2] * k2 + A[6, 3] * k3 + A[6, 4] * k4
                                     + A[ 6, 5] * k5, theta[i - 1], coeff)

        phi[i] = phi[i - 1] + B[1] * k1 + B[2] * k2 + B[3] * k3 + B[4] * k4 + B[5] * k5 + B[6] * k6

        # Solve second differential equation.
        # There is no xi given, so only the y value of the function is altered

        theta[i] = theta[i - 1] + stepSize * phi[i]

    return xi, phi, theta
