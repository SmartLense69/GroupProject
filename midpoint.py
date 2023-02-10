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

def midpoint(h, n, func1, ic1, func2=None, ic2=None):
    """

    :param h:
    :param n:
    :param func1:
    :param ic1:
    :param func2:
    :param ic2:
    :return:
    """

    xivalues = np.arange(0, 1+h, h)
    steps = len(xivalues)
    thetasol = np.arange(0, 1+h, h)
    phisol = np.arange(0, 1+h, h)
    # defining inital conditions
    thetasol[0] = ic2
    phisol[0] = ic1

    for i in range(1, steps):
        # euler to find midpoint
        phisolmid = phisol[i-1] + (h/2)*func1(xivalues[i-1], n, thetasol[i-1], phisol[i-1])
        # using euler to find midpoint
        thetasolmid = thetasol[i - 1] + (h/2) * func2(phisol[i - 1])

        # using midpoint to find next point
        phisol[i] = phisol[i-1] + h*func1(xivalues[i-1], n, thetasolmid, phisolmid)
        # using gradient of midpoint to find next point
        thetasol[i] = thetasol[i-1] + h*func2(phisolmid)

    return xivalues, thetasol


# for l in np.arange(0, 5.5, 0.5):
#     out = midpoint(0.001, l, phigrad, 0, thetagrad, 1)
#     plt.plot(out[0], out[1], label="n=%.1f" %l)
#
#
# plt.legend()
# plt.title('Midpoint Lane-Emden')
# plt.xlabel(r'$\xi$')
# plt.ylabel(r'$\theta$')
# plt.show()
