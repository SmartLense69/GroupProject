import numpy as np
import matplotlib.pyplot as plt

def phigrad(xi, n, theta0, phi0):

    if xi == 0:
        out = 0
    else:
        out = ((-2/xi)*phi0) - theta0**n
    return out


def thetagrad(phi0):
    out = phi0
    return out

def rknorder(stop, h, n, func1, ic1, func2, ic2):
    """

    :param stop: xi value the function runs to
    :param h: step size
    :param n: polytropic n parameter
    :param func1: phigrad
    :param ic1: initial condition for phi
    :param func2: thetagrad
    :param ic2: initial condition for theta
    :return: xi, thetasolution arrays
    """

    xivalues = np.arange(0, stop + h, h)
    steps = len(xivalues)
    thetasol = np.arange(0, stop + h, h)
    phisol = np.arange(0, stop + h, h)
    # defining inital conditions
    thetasol[0] = ic2
    phisol[0] = ic1

    for i in range(1, steps):

        yprime0 = phisol[i-1]
        y0 = thetasol[i-1]
        t0 = xivalues[i-1]

        k1 = func1(t0, n, y0, yprime0)

        yprime1 = yprime0 + k1*h*0.5
        y1 = y0 + 0.5*h*((yprime0+yprime1)/2)

        k2 = func1(t0 + 0.5*h, n, y1, yprime1)

        yprime2 = yprime0 + k2*h*0.5
        y2 = y0 + 0.5*h*((yprime0 + yprime2)/2)

        k3 = func1(t0 + 0.5*h, n, y2, yprime2)

        yprime3 = yprime0 + k3*h
        y3 = y0 + h*((yprime0 + yprime3)/2)

        k4 = func1(t0 + h, n, y3, yprime3)

        phisol[i] = phisol[i-1] + h*(k1 + 2*k2 + 2*k3 + k4)/6
        thetasol[i] = thetasol[i-1] + h*(yprime0 + 2*yprime1 + 2*yprime2 + yprime3)/6

    return xivalues, thetasol

# for k in [0, 1, 5]:
#     output = rknorder(1, 0.001, k, phigrad, 0, thetagrad, 1)
#     plt.plot(output[0], output[1], label='n={0}'.format(k))
#
# plt.legend()
# plt.show()