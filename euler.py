
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


def euler(h, n, func1, ic1, func2=None, ic2=None):
    """
    Function solve the Lane-Emden equation for a given n

    :param h: step size
    :param n: polytropic index
    :param func1: second derivative given as a first derivative with another equation
    :param ic1: initial condition for the second derivative
    :param func2: first derivative given as a constant function
    :param ic2: initial condition for first derivative
    :return: array of xi values and their corresponding theta values
    """

    xivalues = np.arange(0, 1+h, h)
    steps = len(xivalues)
    thetasol = np.arange(0, 1+h, h)
    phisol = np.arange(0, 1+h, h)
    # defining inital conditions
    thetasol[0] = ic2
    phisol[0] = ic1

    for i in range(1, steps):
        if func2 is not None:
            phisol[i] = phisol[i-1] + h*func1(xivalues[i-1], n, thetasol[i-1], phisol[i-1])
        thetasol[i] = thetasol[i-1] + h*func2(phisol[i-1])

    return xivalues, thetasol


# output0 = euler(0.001, 0, phigrad, 0, thetagrad, 1)
# output1 = euler(0.001, 1, phigrad, 0, thetagrad, 1)
# output5 = euler(0.001, 5, phigrad, 0, thetagrad, 1)
# plt.plot(output0[0], output0[1], label='n=0')
# plt.plot(output1[0], output1[1], label='n=1')
# plt.plot(output5[0], output5[1], label='n=5')
# plt.legend()
# plt.title('Euler Lane-Emden')
# plt.xlabel(r'$\xi$')
# plt.ylabel(r'$\theta$')
# plt.show()

def plot(n):
    xiValues, thetaSol = euler(0.001, n, phigrad, 0, thetagrad, 1)
    nLabel = "$n=" + str(n) + "$ with Euler"
    plt.plot(xiValues, thetaSol, label=nLabel)
