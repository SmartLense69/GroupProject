import eulertovdowhile as evw
import numpy as np
import matplotlib.pyplot as plt
import rungekutta4 as rk4
from scipy.interpolate import CubicSpline

CS: CubicSpline = None

def _stateFunction(K, n, P=None, rho=None):
    if P != None:
        return (P / K) ** (n / (n + 1))
    if rho != None:
        return K * (rho ** (1 + 1 / n))


def _massContinuity(r, P, K, n):
    return 4 * np.pi * r ** 2 * _stateFunction(K, n, P=P)


def _hydro(r, P, m, G, c, K, n):
    return - (G * m * _stateFunction(K, n, P=P) / r ** 2)


def _tov(r, P, m, G, c, K, n):
    return - ((G * m * _stateFunction(K, n, P=P)) / r ** 2) \
        * (1 + P / (_stateFunction(K, n, P=P) * c ** 2)) \
        * (1 + (4 * np.pi * r ** 3 * P) / (m * c ** 2)) \
        * (1 - (2 * G * m) / (r * c ** 2)) ** (-1)


def _pressureFunction(rho):
    _x = (1.0088e-2)*np.cbrt(rho/2)
    return 1/(8*(np.pi**2))*(_x * ((2*_x**2)/(3) - 1)
                             * np.sqrt(_x**2 + 1) +
                             np.log(_x + np.sqrt(1 + _x**2)))


def _getRhoSpline(rhoMin=1e+6, rhoMax=1e+9, rhoNum=1000):
    _rho = np.linspace(rhoMin, rhoMax, rhoNum)
    _pressure = _pressureFunction(_rho)

    # Switch x and f of a function f(x),
    # so cs returns x for a certain f(x)
    _cubicSpline = CubicSpline(_pressure, _rho, extrapolate=False)
    return _cubicSpline


def _getRho(P, cubicSpline=None, *args):
    if cubicSpline is not None:
        return cubicSpline(P)
    else:
        cubicSpline = _getRhoSpline(*args)
        return cubicSpline(P)


def _massContinuity2(r, P, K=None, n=None, cs=CS, *args):
    return 4 * np.pi * r ** 2 * _getRho(P, cs)


def _tov2(r, P, m, G, c, K=None, n=None, cs=CS, *args):
    return - ((G * m * cs(P)) / r ** 2) \
        * (1 + P / (cs(P) * c ** 2)) \
        * (1 + (4 * np.pi * r ** 3 * P) / (m * c ** 2)) \
        * (1 - (2 * G * m) / (r * c ** 2)) ** (-1)


def plotMassRadius(rhoMin=1e6, rhoMax=1e9, rhoNum=100):
    _N = 1.5

    plt.rcParams['font.size'] = '14'
    plt.title("Mass radius for various densities")

    _rhoValues = np.linspace(rhoMin, rhoMax, rhoNum)
    _size = _rhoValues.size
    # Index 0 is radius
    # Index 1 is mass

    _massRadiusArray = np.zeros((2, _size))
    cubicSpline = _getRhoSpline(rhoMin, rhoMax, rhoNum)
    CS = cubicSpline

    for i in range(_size):
        _r, _P, _m = rk4.rungekutta4(_massContinuity2, _tov2, 2e10,
                                     1, _rhoValues[i], 6.67e-8, 3e10, 1e13, _N, 1e7, CS)

        _massRadiusArray[0, i] = _r[-1]
        _massRadiusArray[1, i] = _m[-1]

    print(_massRadiusArray)
    plt.plot(_massRadiusArray[0]/100000, _massRadiusArray[1]/(2*1e+33), label="Non-Relativistic\nEquation of state")
    plt.xlabel("Radius in cm")
    plt.ylabel("Mass in g")
    plt.legend()
    plt.grid()
    plt.show()

