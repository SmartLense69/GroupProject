import numpy as np
import matplotlib.pyplot as plt
import euler as eu
import rungekutta4lane as rk4
import rungekuttafehlberg as rukufe

# defining solved analytical solutions
# for n = 0
def analytical0(xi):
    return 1 - ((xi**2)/6)

# for n = 1
def analytical1(xi):
    return np.sin(xi)/xi

# for n = 5
def analytical5(xi):
    return 1/np.sqrt(1+((xi**2)/3))

# defining xi range for the analytical solutions
anaxi = np.arange(1e-12, 1, 0.0001)

h = 0.001

# getting results from the euler scheme
euler0 = eu.euler(5, h, 0, eu.phigrad, 0, eu.thetagrad, 1)
euler1 = eu.euler(5, h, 1, eu.phigrad, 0, eu.thetagrad, 1)
euler5 = eu.euler(5, h, 5, eu.phigrad, 0, eu.thetagrad, 1)

# getting results fron RK4 scheme
rkfour0 = rk4.rungekutta4lane(rk4.phigrad, 2.7, 0, 1, 0, h)
rkfour1 = rk4.rungekutta4lane(rk4.phigrad, 2.7, 0, 1, 1, h)
rkfour5 = rk4.rungekutta4lane(rk4.phigrad, 2.7, 0, 1, 5, h)

# getting results from RKF scheme
diffEq1 = rukufe.de.DifferentialEquation(rukufe.diff1, 0)
diffEq2 = rukufe.de.DifferentialEquation(rukufe.diff2, 1)

rkf0 = rukufe.rkf(diffEq1, diffEq2, h, 3, 0)
rkf1 = rukufe.rkf(diffEq1, diffEq2, h, 3, 1)
rkf5 = rukufe.rkf(diffEq1, diffEq2, h, 3, 5)

# plotting results
plt.figure(figsize=(10, 4), layout='constrained')
plt.suptitle('Lane-Emden Known Analytical Solutions Compared to Scheme Found Solutions')

# for n=0
plt.subplot(131)
plt.plot(anaxi, analytical0(anaxi), label='Analytical', linestyle='dashed')
plt.plot(euler0[0], euler0[1], label='Euler')
plt.plot(rkfour0[0], rkfour0[1], label='RK4')
rukufe.plot(0)
# plt.plot(rkf0[0], rkf0[1], label='RKF')
plt.xlim([0, 1])
plt.ylim([0, 1.1])
plt.title('n = 0')
plt.legend()
plt.xlabel(r'$\xi$')
plt.ylabel(r'$\theta$')

plt.subplot(132)
plt.plot(anaxi, analytical1(anaxi), label='Analytical', linestyle='dashed')
plt.plot(euler1[0], euler1[1], label='Euler')
plt.plot(rkfour1[0], rkfour1[1], label='RK4')
rukufe.plot(1)
# plt.plot(rkf1[0], rkf1[1], label='RKF')
plt.xlim([0, 1])
plt.ylim([0, 1.1])
plt.title('n = 1')
plt.legend()
plt.xlabel(r'$\xi$')
plt.ylabel(r'$\theta$')

plt.subplot(133)
plt.plot(anaxi, analytical5(anaxi), label='Analytical', linestyle='dashed')
plt.plot(euler5[0], euler5[1], label='Euler')
plt.plot(rkfour5[0], rkfour5[1], label='RK4')
rukufe.plot(5)
# plt.plot(rkf5[0], rkf5[1], label='RKF')
plt.xlim([0, 1])
plt.ylim([0, 1.1])
plt.title('n = 5')
plt.legend()
plt.xlabel(r'$\xi$')
plt.ylabel(r'$\theta$')

plt.show()

# def rkf(diffEq1: de.DifferentialEquation, diffEq2: de.DifferentialEquation, stepSize, maxTime=1, coeff=1):
#     """
#     Solves a differential equation system.
#     This function only applies for Lane-Emden currently.
#     :param diffEq1: First decoupled diff. equation (-theta^4+2phi/xi)
#     :param diffEq2: Second decoupled diff. equation (phi)
#     :param stepSize: Size of the step, often mentioned in the scripts as h
#     :param maxTime: The max time the differential equation should be solved to.
#     :param coeff: Additional parameter for Lane-Emden
#     :return: xi, phi and theta of the set of differential equations
#     """
#
#     xi = np.arange(0, maxTime + stepSize, stepSize)
#     n = np.shape(xi)[0]
#     phi = np.zeros(n)
#     theta = np.zeros(n)
#     phi[0] = diffEq1.iniC
#     theta[0] = diffEq2.iniC
#
#     for i in range(1, n):
#
#         # Solve the first differential equation
#         # Keep theta constant
#
#         k1 = stepSize * diffEq1.func(xi[i], phi[i - 1], theta[i - 1], coeff)
#         k2 = stepSize * diffEq1.func(xi[i] + C[2] * stepSize,
#                                      phi[i - 1] + A[2, 1] * k1, theta[i - 1], coeff)
#         k3 = stepSize * diffEq1.func(xi[i] + C[3] * stepSize,
#                                      phi[i - 1] + A[3, 1] * k1 + A[3, 2] * k2, theta[i - 1], coeff)
#         k4 = stepSize * diffEq1.func(xi[i] + C[4] * stepSize,
#                                      phi[i - 1] + A[4, 1] * k1 + A[4, 2] * k2 + A[4, 3] * k3, theta[i - 1], coeff)
#         k5 = stepSize * diffEq1.func(xi[i] + C[5] * stepSize,
#                                      phi[i - 1] + A[5, 1] * k1 + A[5, 2] * k2 + A[5, 3] * k3 + A[5, 4] * k4,
#                                      theta[i - 1], coeff)
#         k6 = stepSize * diffEq1.func(xi[i] + C[6] * stepSize,
#                                      phi[i - 1] + A[6, 1] * k1 + A[6, 2] * k2 + A[6, 3] * k3 + A[6, 4] * k4
#                                      + A[ 6, 5] * k5, theta[i - 1], coeff)
#
#         phi[i] = phi[i - 1] + B[1] * k1 + B[2] * k2 + B[3] * k3 + B[4] * k4 + B[5] * k5 + B[6] * k6
#
#         # Solve second differential equation.
#         # There is no xi given, so only the y value of the function is altered
#
#         theta[i] = theta[i - 1] + stepSize * phi[i]
#
#     return xi, phi, theta