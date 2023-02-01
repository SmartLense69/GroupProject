import numpy as np
import matplotlib.pyplot as plt
import diffequation as de
import euler_CR as eCR

# np.seterr(all='raise')


def _phigrad(xi, theta0, phi0):
    if xi == 0:
        return xi
    else:
        return ((-2/xi)*phi0) - theta0


def _thetagrad(phi0):
    return phi0


def euler(stop, h, n, func1, ic1, func2=None, ic2=None):
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

    xivalues = np.arange(0, stop+h, h)
    steps = len(xivalues)
    thetasol = np.arange(0, stop+h, h)
    phisol = np.arange(0, stop+h, h)
    # defining inital conditions
    thetasol[0] = ic2
    phisol[0] = ic1

    for i in range(1, steps):
        if func2 is not None:
            phisol[i] = phisol[i-1] + h*func1(xivalues[i-1], n, thetasol[i-1], phisol[i-1])
        thetasol[i] = thetasol[i-1] + h*func2(phisol[i-1])

    return xivalues, thetasol


# output0 = euler(5, 0.001, 0, phigrad, 0, thetagrad, 1)
# output1 = euler(5, 0.001, 1, phigrad, 0, thetagrad, 1)
# output5 = euler(5, 0.001, 5, phigrad, 0, thetagrad, 1)
# plt.plot(output0[0], output0[1], label='n=0')
# plt.plot(output1[0], output1[1], label='n=1')
# plt.plot(output5[0], output5[1], label='n=5')

# goto = 35
# for n in np.arange(0, 5.5, 0.5):
#     xiValues, thetaSol = euler(goto, 0.001, n, phigrad, 0, thetagrad, 1)
#     nLabel = "$n=" + str(n) + "$ with Euler"
#     plt.plot(xiValues, thetaSol, label=nLabel)
#
# plt.legend()
# plt.hlines(0, 0, goto, color='black', linestyles='--')
# plt.title('Euler Lane-Emden')
# plt.xlabel(r'$\xi$')
# plt.ylabel(r'$\theta$')
# plt.ylim([-1, 1.2])
# plt.show()

def plot(n):
    diff1 = de.DifferentialEquation(_phigrad, 0)
    diff2 = de.DifferentialEquation(_thetagrad, 1)

    _xiValues, _phiValues, _thetaValues = eCR.euler2Diff(diff1, diff2, 0.001, 5)

    # xiValues, thetaSol = euler(5, 0.001, n, _phigrad, 0, _thetagrad, 1)
    nLabel = "$n=" + str(n) + "$ with Euler"
    plt.plot(_xiValues, _thetaValues, label=nLabel)