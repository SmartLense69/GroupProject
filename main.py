import DiffEquation as DE
import RKF
import matplotlib.pyplot as plt

N = 1


def diff1(xi, phi, theta):
    return -theta**N - (2/xi) * phi


def diff2(phi):
    return phi


if __name__ == '__main__':
    diffEq1 = DE.DifferentialEquation(diff1, 0)
    diffEq2 = DE.DifferentialEquation(diff2, 1)
    xi, phi, theta = RKF.rkf(diffEq1, diffEq2, 0.01, 1)

    plt.title("Solution to Lane-Emden")
    plt.plot(xi, theta)
    plt.show()