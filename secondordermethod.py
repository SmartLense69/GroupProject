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

def secondorder(stop, h, n, func1, ic1, func2, ic2):
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

        j1 = h*func1(xivalues[i-1], n, thetasol[i-1], phisol[i-1])
        k1 = h*func2(phisol[i-1])

        j2 = h*func1(xivalues[i-1] + 0.5*h, n, thetasol[i-1] + 0.5*k1, phisol[i-1] + 0.5*j1)
        k2 = h*func2(phisol[i-1] + 0.5*j1)

        phisol[i] = phisol[i-1] + j2
        thetasol[i] = thetasol[i-1] + k2

        # phisol[i] = phisol[i - 1] + h * func1(xivalues[i - 1] + 0.5*h, n, thetasol[i - 1], phisol[i - 1]+ 0.5*h*func1(xivalues[i-1], n, thetasol[i-1], phisol[i-1]))
        # thetasol[i] = thetasol[i - 1] + h * func2(phisol[i - 1])

    return xivalues, thetasol

for k in [0, 1, 5]:
    output = secondorder(5, 0.001, k, phigrad, 0, thetagrad, 1)
    plt.plot(output[0], output[1], label='n=0')

plt.show()