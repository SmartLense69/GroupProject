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

def backwardeuler(stop, h, n, func1, ic1, func2, ic2):
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
        # phisol[i] = phisol[i - 1] + h * func1(xivalues[i], n, thetasol[i], phisol[i])

        phisol[i] - phisol[i-1] - h*func1(xivalue[i], n, thetalsol[i-1], phisol[i]) = 0

        # thetasol[i] = thetasol[i - 1] + h * func2(phisol[i])


    return xivalues, thetasol

for k in [0, 1, 5]:
    output = backwardeuler(1, 0.001, k, phigrad, 0, thetagrad, 1)
    plt.plot(output[0], output[1], label='n={0}'.format(k))

plt.legend()
plt.show()