import diffequation as de
import rungekutta5th as rk5
import matplotlib.pyplot as plt
import numpy as np

N = 0


def diff1(inXi, inPhi, inTheta):
    return -inTheta**N - (2/inXi) * inPhi


def diff2(inPhi):
    return inPhi


if __name__ == '__main__':

    # Define two differential equations, which will
    # represent the two decoupled Lane-Emden eqs.

    diffEq1 = de.DifferentialEquation(diff1, 0)
    diffEq2 = de.DifferentialEquation(diff2, 1)

    # Solve with RKF

    xi0, phi0, theta0 = rk5.rkf(diffEq1, diffEq2, 0.01, 1)

    N = 1
    xi1, phi1, theta1 = rk5.rkf(diffEq1, diffEq2, 0.01, 1)

    N = 5
    xi5, phi5, theta5 = rk5.rkf(diffEq1, diffEq2, 0.01, 1)

    # Plot
    plt.title("Solution to Lane-Emden for $n=0,1,5$")
    plt.plot(xi0, theta0, label="$n=0$")
    plt.plot(xi1, theta1, label="$n=1$")
    plt.plot(xi5, theta5, label="$n=5$")
    plt.grid()
    plt.legend()
    plt.xlim(0, 1)
    plt.xlabel("Non-dimensional radius from 0 to 1")
    plt.ylabel("Non-dimensional density from 0 to 1")
    plt.show()

    # Iterate through multiple N

    NArray = np.arange(start=1, stop=5.5, step=0.5)
    plt. title("Solution to Lane-Emden for $1<n<5$")

    for n in NArray:
        N = n
        xi, phi, theta = rk5.rkf(diffEq1, diffEq2, 0.01, 1)
        nLabel = "$n=" + str(N) + "$"
        plt.plot(xi0, theta0, label=nLabel)

    plt.grid()
    plt.legend()
    plt.xlabel("Non-dimensional radius from 0 to 1")
    plt.ylabel("Non-dimensional density from 0 to 1")
    plt.show()
