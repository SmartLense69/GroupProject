import DiffEquation as DE
import RKF
import matplotlib.pyplot as plt

N = 1.5


def diff1(xi, phi, theta):
    return -theta**N - (2/xi) * phi


def diff2(phi):
    return phi


if __name__ == '__main__':
    diffEq1 = DE.DifferentialEquation(diff1, 1)
    diffEq2 = DE.DifferentialEquation(diff2, 0)
    xi, theta, phi = RKF.rkf(diffEq1, diffEq2, 0.1, 1)

    plt.title("Solution to Lane-Emden")
    plt.plot(xi, theta)
    plt.show()