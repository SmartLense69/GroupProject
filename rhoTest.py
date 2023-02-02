import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline


# Global variables

_M0 = 1
_G = 6.67e-8
_C = 3e10
_K = 1e13
_N = 1.5
_R = 2e10
_PLOTPRHO = True

cubicSplineData = np.loadtxt("EoS_Spline_Data.csv", delimiter='\t').transpose()
getRho = CubicSpline(cubicSplineData[1], cubicSplineData[0], extrapolate=True)
getP = CubicSpline(cubicSplineData[0], cubicSplineData[1], extrapolate=True)

if _PLOTPRHO:
    _p = np.linspace(1e15, 1e19, 1000)
    plt.plot(_p, getRho(_p))
    plt.show()


def _M(r, P):
    return 4 * np.pi * r ** 2 * getRho(P)


def _TOV(r, P, m):
    return - ((_G * m * getRho(P)) / r ** 2) \
        * (1 + P / (getRho(P) * _C ** 2)) \
        * (1 + (4 * np.pi * r ** 3 * P) / (m * _C ** 2)) \
        * (1 - (2 * _G * m) / (r * _C ** 2)) ** (-1)


def _RK4(masscont, starfunc, statefunc, h, rho):

    # constants
    c = np.array([None, 0, 1/2, 1/2, 1])
    a = np.array([[None, None, None, None],
                 [None, None, None, None],
                 [None, 1/2, None, None],
                 [None, 0, 1/2, None],
                 [None, 0, 0, 1]])
    b = np.array([None, 1/6, 1/3, 1/3, 1/6])

    # initialise the arrays to be used

    r = np.arange(1, _R + h, h)
    # number of time steps
    steps = np.shape(r)[0]

    m = np.zeros(steps)
    P = np.zeros(steps)

    # set the initial conditions
    m[0] = _M0
    P[0] = statefunc(rho)

    # Loop over time
    for i in range(1, steps):

        # find m using RK4 method
        # calculate the values k1, j1 through k4, j4
        j1 = h * masscont(r[i - 1] + c[1] * h, np.abs(P[i - 1]))
        k1 = h * starfunc(r[i - 1] + c[1] * h, np.abs(P[i - 1]), m[i - 1])

        j2 = h * masscont(r[i - 1] + c[2] * h, np.abs(P[i - 1]) + a[2, 1] * k1)
        k2 = h * starfunc(r[i - 1] + c[2] * h, np.abs(P[i - 1]) + a[2, 1] * k1, m[i - 1] + a[2, 1] * j1)

        j3 = h * masscont(r[i - 1] + c[3] * h, np.abs(P[i - 1]) + a[3, 1] * k1 + a[3, 2] * k2)
        k3 = h * starfunc(r[i - 1] + c[3] * h, np.abs(P[i - 1]) + a[3, 1] * k1 + a[3, 2] * k2,
                          m[i - 1] + a[3, 1] * j1 + a[3, 2] * j2)

        j4 = h * masscont(r[i - 1] + c[4] * h, np.abs(P[i - 1]) + a[4, 1] * k1 + a[4, 2] * k2 + a[4, 3] * k3)
        k4 = h * starfunc(r[i - 1] + c[4] * h, np.abs(P[i - 1]) + a[4, 1] * k1 + a[4, 2] * k2 + a[4, 3] * k3,
                          m[i - 1] + a[4, 1] * j1 + a[4, 2] * j2 + a[4, 3] * j3)

        # find next value for m
        m[i] = m[i - 1] + b[1] * j1 + b[2] * j2 + b[3] * j3 + b[4] * j4

        # find next value for P
        P[i] = np.abs(P[i - 1]) + b[1] * k1 + b[2] * k2 + b[3] * k3 + b[4] * k4

        if P[i]/P[0] < 1e-5:
            P = P[:i]
            r = r[:i]
            m = m[:i]
            break

    return r, P, m


def runRK4():
    for rho in np.geomspace(1e6, 1e12, 100):
        r, P, m = _RK4(_M, _TOV, getP, 1e5, rho)
        plt.scatter(r[-1], m[-1], color='b', marker='v')
    plt.grid()
    plt.show()

