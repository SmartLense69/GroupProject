import eulertovdowhile as evw
import numpy as np
import matplotlib.pyplot as plt
import rungekutta4 as rk4


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


def _pressureFunction(rho, K, n, P=None):
    _x = (1.0088e-2)*np.cbrt(rho/2)
    return 1/(8*(np.pi**2))*(_x * ((2*_x**2)/(3) - 1)
                             * np.sqrt(_x**2 + 1) +
                             np.ln(_x + np.sqrt(1 + _x**2)))


def plotMassRadius(rhoMin=1e6, rhoMax=1e9, rhoNum=1000):
    _N = 1.5

    plt.rcParams['font.size'] = '14'
    plt.title("Mass radius for various densities")

    _rhoValues = np.linspace(rhoMin, rhoMax, rhoNum)
    _size = _rhoValues.size
    # Index 0 is radius
    # Index 1 is mass

    _massRadiusArray = np.zeros((2, _size))
    for i in range(_size):
        _r, _P, _m = rk4.rungekutta4(_massContinuity, _tov, 2e10,
                                     1, _rhoValues[i], 6.67e-8, 3e10, 1e13, _N, 1e5)

        _massRadiusArray[0, i] = _r[-1]
        _massRadiusArray[1, i] = _m[-1]

    print(_massRadiusArray)
    plt.plot(_massRadiusArray[0]/100000, _massRadiusArray[1]/(2*1e+33), label="Non-Relativistic\nEquation of state")
    plt.xlabel("Radius in cm")
    plt.ylabel("Mass in g")
    plt.legend()
    plt.grid()
    plt.show()
