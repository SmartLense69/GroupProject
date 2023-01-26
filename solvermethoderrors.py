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

# # plotting results
# plt.figure(figsize=(10, 4), layout='constrained')
# plt.suptitle('Lane-Emden Known Analytical Solutions Compared to Scheme Found Solutions')
#
# # for n=0
# plt.subplot(131)
# plt.plot(anaxi, analytical0(anaxi), label='Analytical', linestyle='dashed')
# plt.plot(euler0[0], euler0[1], label='Euler')
# plt.plot(rkfour0[0], rkfour0[1], label='RK4')
# rukufe.plot(0)
# # plt.plot(rkf0[0], rkf0[1], label='RKF')
# plt.xlim([0, 1])
# plt.ylim([0, 1.1])
# plt.title('n = 0')
# plt.legend()
# plt.xlabel(r'$\xi$')
# plt.ylabel(r'$\theta$')
#
# plt.subplot(132)
# plt.plot(anaxi, analytical1(anaxi), label='Analytical', linestyle='dashed')
# plt.plot(euler1[0], euler1[1], label='Euler')
# plt.plot(rkfour1[0], rkfour1[1], label='RK4')
# rukufe.plot(1)
# # plt.plot(rkf1[0], rkf1[1], label='RKF')
# plt.xlim([0, 1])
# plt.ylim([0, 1.1])
# plt.title('n = 1')
# plt.legend()
# plt.xlabel(r'$\xi$')
# plt.ylabel(r'$\theta$')
#
# plt.subplot(133)
# plt.plot(anaxi, analytical5(anaxi), label='Analytical', linestyle='dashed')
# plt.plot(euler5[0], euler5[1], label='Euler')
# plt.plot(rkfour5[0], rkfour5[1], label='RK4')
# rukufe.plot(5)
# # plt.plot(rkf5[0], rkf5[1], label='RKF')
# plt.xlim([0, 1])
# plt.ylim([0, 1.1])
# plt.title('n = 5')
# plt.legend()
# plt.xlabel(r'$\xi$')
# plt.ylabel(r'$\theta$')
#
# plt.show()

# errors with varying h


# n = 0
testh = [1, 1e-1, 5e-1, 1e-2, 5e-2, 1e-3, 5e-3, 1e-4, 5e-4, 1e-5, 5e-5, 1e-6]
# testh = np.arange(1e-6, 1e-1, 5e-5)
eulererror0 = np.zeros(len(testh))
rkfourerror0 = np.zeros(len(testh))
count = np.arange(0, len(testh), 1)
for i in count:
    h = testh[i]
    erroreuler0 = eu.euler(1, h, 0, eu.phigrad, 0, eu.thetagrad, 1)
    errorrk40 = rk4.rungekutta4lane(rk4.phigrad, 1, 0, 1, 0, h)
    abseuler0 = np.abs(analytical0(erroreuler0[0][len(erroreuler0[0])-1]) - erroreuler0[1][len(erroreuler0[1])-1])
    absrkfour0 = np.abs(analytical0(errorrk40[0][len(errorrk40[0])-1]) - errorrk40[1][len(errorrk40[1])-1])
    eulererror0[i] = abseuler0
    rkfourerror0[i] = absrkfour0

plt.scatter(testh, eulererror0, label='Euler')
plt.scatter(testh, rkfourerror0, label='RK4')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('h')
plt.ylabel('Absolute Error')
plt.title('Running to 1')
plt.show()