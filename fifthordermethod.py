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

def fifthorder(stop, h, n, func1, ic1, func2, ic2):
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

        j1 = h*func1(xivalues[i-1] + 0*h, n, thetasol[i-1], phisol[i-1])
        k1 = h*func2(phisol[i-1])

        j2 = h*func1(xivalues[i-1] + h/4, n, thetasol[i-1] + k1/4 , phisol[i-1] + j1/4)
        k2 = h*func2(phisol[i-1] + j1/4)

        j3 = h*func1(xivalues[i-1] + 3*h/8, n, thetasol[i-1] + 3*k1/32 + 9*k2/32, phisol[i-1] + 3*j1/32 + 9*j2/32)
        k3 = h*func2(phisol[i-1] + 3*j1/32 + 9*j2/32)

        j4 = h*func1(xivalues[i-1] + 12*h/13, n, thetasol[i-1] + 1932*k1/2197 - 7200*k2/2197 + 7296*k3/2197, phisol[i-1] + 1932*j1/2197 - 7200*j2/2197 + 7296*j3/2197)
        k4 = h*func2(phisol[i-1] + 1932*j1/2197 - 7200*j2/2197 + 7296*j3/2197)

        j5 = h*func1(xivalues[i-1] + h, n, thetasol[i-1] + 439*k1/216 - 8*k2 + 3680*k3/513 - 845*k4/4104, phisol[i-1] + 439*j1/216 - 8*j2 + 3680*j3/513 - 845*j4/4104)
        k5 = h*func2(phisol[i-1] + 439*j1/216 - 8*j2 + 3680*j3/513 - 845*j4/4104)

        j6 = h*func1(xivalues[i-1] + 0.5*h, n, thetasol[i-1] - 8*k1/27 + 2*k2 - 3544*k3/2565 + 1859*k4/4104 - 11*k5/40, phisol[i-1] - 8*j1/27 + 2*j2 - 3544*j3/2565 + 1859*j4/4104 - 11*j5/40)
        k6 = h*func2(phisol[i-1] - 8*j1/27 + 2*j2 - 3544*j3/2565 + 1859*j4/4104 - 11*j5/40)

        phisol[i] = phisol[i-1] + 16*j1/135 + 0*j2 + 6656*j3/12825 + 28561*j4/56430 - 9*j5/50 + 2*j6/55
        thetasol[i] = thetasol[i-1] + 16*k1/135 + 0*k2 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55


    return xivalues, thetasol

for k in [0, 1, 5]:
    output = fifthorder(1, 0.001, k, phigrad, 0, thetagrad, 1)
    plt.plot(output[0], output[1], label='n={0}'.format(k))

plt.legend()
plt.show()