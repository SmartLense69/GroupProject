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


def _pressureFunction(K=None, n=None, rho=None):
    _x = (1.0088e-2)*np.cbrt(rho/2)
    _phi = 1/(8*(np.pi**2))*(_x * ((2*_x**2)/(3) - 1)
                             * np.sqrt(_x**2 + 1) +
                             np.log(_x + np.sqrt(1 + _x**2)))
    _P = 1.42180e25 * _phi
    return _P


def generateSplineData(rhoMinExp=2, rhoMaxExp=15, rhoAdaptStepNum=100):
    csvArray = np.zeros(shape=((rhoMaxExp - rhoMinExp) * 100, 2))
    i = 0
    for n in range(rhoMinExp, rhoMaxExp):
        _rho = np.linspace(1 * (10 ** (n - 2)), 1 * (10 ** n), rhoAdaptStepNum)
        _pressure = _pressureFunction(rho=_rho)
        for rho, p in zip(_rho, _pressure):
            if i >= csvArray.shape[0]:
                i = 0
                break
            csvArray[i, 0] = rho
            csvArray[i, 1] = p
            print(i)
            i += 1
    np.savetxt("EoS_Spline_Data.csv", csvArray, delimiter='\t')

def _getRhoSpline(rhoMin=1e+6, rhoMax=1e+9, rhoNum=100):
    _rho = np.linspace(rhoMin, rhoMax, rhoNum)
    _pressure = _pressureFunction(rho=_rho)

    # Switch x and f of a function f(x),
    # so cs returns x for a certain f(x)
    _cubicSpline = CubicSpline(_pressure, _rho, extrapolate=True)
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


def plotMassRadius(rhoMin=1e6, rhoMax=5e10, rhoNum=1000):
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
        _r, _P, _m = rk4.rungekutta4(_massContinuity2, _tov2, _pressureFunction, 2e10,
                                     1, _rhoValues[i], 6.67e-8, 3e10, 1e13, _N, 1e6, CS)

        _massRadiusArray[0, i-1] = _r[-1]
        _massRadiusArray[1, i-1] = _m[-1]

    print(_massRadiusArray)
    plt.plot(_massRadiusArray[0, :-1]/100000, _massRadiusArray[1, :-1]/(2*1e+33), label="Non-Relativistic\nEquation of state")
    plt.xlabel("Radius in cm")
    plt.ylabel("Mass in g")
    plt.legend()
    plt.grid()
    plt.show()

