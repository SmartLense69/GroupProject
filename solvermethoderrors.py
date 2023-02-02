import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
import euler as eu
import rungekutta4lane as rk4l
import rungekutta4 as rk4
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


# # absolute errors with varying h
# # range of step sizes to test
# testh = [1, 1e-1, 5e-1, 1e-2, 5e-2, 1e-3]#, 5e-3, 1e-4, 5e-4] #, 1e-5, 5e-5, 1e-6]
# # testh = np.arange(1e-6, 1e-1, 5e-5)
#
# # defining solution arrays for the absolute error values
# eulererror0 = np.zeros(len(testh))
# rkfourerror0 = np.zeros(len(testh))
# eulererror1 = np.zeros(len(testh))
# rkfourerror1 = np.zeros(len(testh))
# eulererror5 = np.zeros(len(testh))
# rkfourerror5 = np.zeros(len(testh))
# rkferror0 = np.zeros(len(testh))
# rkferror1 = np.zeros(len(testh))
# rkferror5 = np.zeros(len(testh))
# count = np.arange(0, len(testh), 1)
#
# # loop for calculating errors for each step size value
# for i in count:
#     # setting h for the iteration of the loop
#     h = testh[i]
#
#     # n=0
#     # finding the euler value at xi = 1
#     erroreuler0 = eu.euler(1, h, 0, eu.phigrad, 0, eu.thetagrad, 1)
#     # finding the rk4 value at xi = 1
#     errorrk40 = rk4.rungekutta4lane(rk4.phigrad, 1, 0, 1, 0, h)
#     # finding the rkf value at xi = 1
#     errorrkf0 = rukufe.rkf(rukufe.de.DifferentialEquation(rukufe.diff1, 0), rukufe.de.DifferentialEquation(rukufe.diff2, 1), h, 1, 0)
#     # calculating the absolute error in euler at xi = 1 by subtracting from the analytical solution
#     abseuler0 = np.abs(analytical0(erroreuler0[0][len(erroreuler0[0])-1]) - erroreuler0[1][len(erroreuler0[1])-1])
#     # calculating the absolute error in rk4 at xi = 1 by subtracting from the analytical solution
#     absrkfour0 = np.abs(analytical0(errorrk40[0][len(errorrk40[0])-1]) - errorrk40[1][len(errorrk40[1])-1])
#     # calculating the absolute error in rkf at xi = 1 by subtracting from the analytical solution
#     absrkf0 = np.abs(analytical0(errorrkf0[0][len(errorrkf0[0])-1]) - errorrkf0[2][len(errorrkf0[2])-1])
#     # adding results to the solution arrays
#     eulererror0[i] = abseuler0
#     rkfourerror0[i] = absrkfour0
#     rkferror0[i] = absrkf0
#
#     # n=1, repeating same steps as the n=0 case
#     erroreuler1 = eu.euler(1, h, 1, eu.phigrad, 0, eu.thetagrad, 1)
#     errorrk41 = rk4.rungekutta4lane(rk4.phigrad, 1, 0, 1, 1, h)
#     errorrkf1 = rukufe.rkf(rukufe.de.DifferentialEquation(rukufe.diff1, 0),
#                            rukufe.de.DifferentialEquation(rukufe.diff2, 1), h, 1, 1)
#     abseuler1 = np.abs(analytical1(erroreuler1[0][len(erroreuler1[0]) - 1]) - erroreuler1[1][len(erroreuler1[1]) - 1])
#     absrkfour1 = np.abs(analytical1(errorrk41[0][len(errorrk41[0]) - 1]) - errorrk41[1][len(errorrk41[1]) - 1])
#     absrkf1 = np.abs(analytical1(errorrkf1[0][len(errorrkf1[0]) - 1]) - errorrkf1[2][len(errorrkf1[2]) - 1])
#     eulererror1[i] = abseuler1
#     rkfourerror1[i] = absrkfour1
#     rkferror1[i] = absrkf1
#
#     # n=5, repeating the same steps as the n=0 case
#     erroreuler5 = eu.euler(1, h, 5, eu.phigrad, 0, eu.thetagrad, 1)
#     errorrk45 = rk4.rungekutta4lane(rk4.phigrad, 1, 0, 1, 5, h)
#     errorrkf5 = rukufe.rkf(rukufe.de.DifferentialEquation(rukufe.diff1, 0),
#                            rukufe.de.DifferentialEquation(rukufe.diff2, 1), h, 1, 5)
#     abseuler5 = np.abs(analytical5(erroreuler5[0][len(erroreuler5[0])-1]) - erroreuler5[1][len(erroreuler5[1])-1])
#     absrkfour5 = np.abs(analytical5(errorrk45[0][len(errorrk45[0]) - 1]) - errorrk45[1][len(errorrk45[1]) - 1])
#     absrkf5 = np.abs(analytical5(errorrkf0[0][len(errorrkf0[0]) - 1]) - errorrkf5[2][len(errorrkf5[2]) - 1])
#     eulererror5[i] = abseuler5
#     rkfourerror5[i] = absrkfour5
#     rkferror5[i] = absrkf5
#
# # fitting the error plots to dependence on h
# # fitting function
# def fitlinear(x, a, b):
#     return a*x + b
#
# # defining range of values to plot fit function over, sufficient to be smooth line
# fith_smooth = np.linspace(min(testh), max(testh), 1000)
#
# # euler fits
# # n = 0
# euler0_params, params_covariance = optimize.curve_fit(fitlinear, np.log10(testh), np.log10(eulererror0), p0=[1, 1])
# eulererrorfit0 = fitlinear(np.log10(fith_smooth),euler0_params[0], euler0_params[1])
# print("Error \u221d h^{0:.3g}".format(euler0_params[0]))
# # n = 1
# euler1_params, params_covariance = optimize.curve_fit(fitlinear, np.log10(testh), np.log10(eulererror1), p0=[1, 1])
# eulererrorfit1 = fitlinear(np.log10(fith_smooth),euler1_params[0], euler1_params[1])
# print("Error \u221d h^{0:.3g}".format(euler1_params[0]))
# # n = 5
# euler5_params, params_covariance = optimize.curve_fit(fitlinear, np.log10(testh), np.log10(eulererror5), p0=[1, 1])
# eulererrorfit5 = fitlinear(np.log10(fith_smooth),euler5_params[0], euler5_params[1])
# print("Error \u221d h^{0:.3g}".format(euler5_params[0]))
#
# # rk4 fits
# # n = 0
# rk40_params, params_covariance = optimize.curve_fit(fitlinear, np.log10(testh), np.log10(rkfourerror0), p0=[1, 1])
# rk4errorfit0 = fitlinear(np.log10(fith_smooth),rk40_params[0], rk40_params[1])
# print("Error \u221d h^{0:.3g}".format(rk40_params[0]))
# # n = 1
# rk41_params, params_covariance = optimize.curve_fit(fitlinear, np.log10(testh), np.log10(rkfourerror1), p0=[1, 1])
# rk4errorfit1 = fitlinear(np.log10(fith_smooth),rk41_params[0], rk41_params[1])
# print("Error \u221d h^{0:.3g}".format(rk41_params[0]))
# # n = 5
# rk45_params, params_covariance = optimize.curve_fit(fitlinear, np.log10(testh), np.log10(rkfourerror5), p0=[1, 1])
# rk4errorfit5 = fitlinear(np.log10(fith_smooth),rk45_params[0], rk45_params[1])
# print("Error \u221d h^{0:.3g}".format(rk45_params[0]))
#
# # rkf fits
# # n = 0
# rkf0_params, params_covariance = optimize.curve_fit(fitlinear, np.log10(testh), np.log10(rkferror0), p0=[1, 1])
# rkferrorfit0 = fitlinear(np.log10(fith_smooth),rkf0_params[0], rkf0_params[1])
# print("Error \u221d h^{0:.3g}".format(rkf0_params[0]))
# # n = 1
# rkf1_params, params_covariance = optimize.curve_fit(fitlinear, np.log10(testh), np.log10(rkferror1), p0=[1, 1])
# rkferrorfit1 = fitlinear(np.log10(fith_smooth),rkf1_params[0], rkf1_params[1])
# print("Error \u221d h^{0:.3g}".format(rkf1_params[0]))
# # n = 5
# rkf5_params, params_covariance = optimize.curve_fit(fitlinear, np.log10(testh), np.log10(rkferror5), p0=[1, 1])
# rkferrorfit5 = fitlinear(np.log10(fith_smooth),rkf5_params[0], rkf5_params[1])
# print("Error \u221d h^{0:.3g}".format(rkf5_params[0]))
#
# # creating plot of all the results
# plt.figure(figsize=(12, 4), layout='constrained')
# plt.suptitle('Error Plots Running to 1')
#
# # plotting the n=0 errors against step size for euler and rk4
# plt.subplot(131)
# plt.scatter(testh, eulererror0, label='Euler')
# plt.plot(fith_smooth, 10**eulererrorfit0, label='Euler Fit')
# plt.scatter(testh, rkfourerror0, label='RK4')
# plt.plot(fith_smooth, 10**rk4errorfit0, label='RK4 Fit')
# plt.scatter(testh, rkferror0, label='RKF')
# plt.plot(fith_smooth, 10**rkferrorfit0, label='RKf Fit')
# plt.legend()
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('h')
# plt.ylabel('Absolute Error')
# plt.title('n = 0')
#
# # plotting the n=1 errors against step size for euler and rk4
# plt.subplot(132)
# plt.scatter(testh, eulererror1, label='Euler')
# plt.plot(fith_smooth, 10**eulererrorfit1, label='Euler Fit')
# plt.scatter(testh, rkfourerror1, label='RK4')
# plt.plot(fith_smooth, 10**rk4errorfit1, label='RK4 Fit')
# plt.scatter(testh, rkferror1, label='RKF')
# plt.plot(fith_smooth, 10**rkferrorfit1, label='RKf Fit')
# plt.legend()
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('h')
# plt.ylabel('Absolute Error')
# plt.title('n = 1')
#
# # plotting the n=5 errors against step size for euler and rk4
# plt.subplot(133)
# plt.scatter(testh, eulererror5, label='Euler')
# plt.plot(fith_smooth, 10**eulererrorfit5, label='Euler Fit')
# plt.scatter(testh, rkfourerror5, label='RK4')
# plt.plot(fith_smooth, 10**rk4errorfit5, label='RK4 Fit')
# plt.scatter(testh, rkferror5, label='RKF')
# plt.plot(fith_smooth, 10**rkferrorfit5, label='RKf Fit')
# plt.legend()
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('h')
# plt.ylabel('Absolute Error')
# plt.title('n = 5')
#
# plt.show()

# rewritting with loops

testh = [1, 1e-1, 5e-1, 1e-2, 5e-2, 1e-3, 5e-3, 1e-4, 5e-4] #, 1e-5, 5e-5, 1e-6]

# fitting function
def fitlinear(x, a, b):
    return a*x + b

# defining range of values to plot fit function over, sufficient to be smooth line
fith_smooth = np.linspace(min(testh), max(testh), 1000)

plt.figure(figsize=(12, 4), layout='constrained')
plt.suptitle('Error Plots Running to 1')

js = [0, 1, 5]
indices = [1, 2, 3]
analyt = [analytical0, analytical1, analytical5]

for j, index, ana in zip(js, indices, analyt):

    # euler
    eulererrors = np.zeros(len(testh))
    rk4errors = np.zeros(len(testh))

    for i in range(0, len(testh)):
        h = testh[i]

        eulerto1 = eu.euler(3, h, j, eu.phigrad, 0, eu.thetagrad, 1)
        abserroreulerat1 = np.abs(ana(eulerto1[0][len(eulerto1[0]) - 1]) - eulerto1[1][len(eulerto1[1]) - 1])
        eulererrors[i] = abserroreulerat1

        #rk4to1 = rk4.rungekutta4lane(rk4.phigrad, 3, 0, 1, j, h)
        rk4to1 = rk4.rungekutta4(rk4.phigrad, rk4.thetagrad, rk4.statefunc, 3, 0, 1, j, h)
        abserrorrk4at1 = np.abs(ana(rk4to1[0][-1]) - rk4to1[2][-1])
        rk4errors[i] = abserrorrk4at1

    euler_params, params_covariance = optimize.curve_fit(fitlinear, np.log10(testh), np.log10(eulererrors), p0=[1, 1])
    eulererrorfit = fitlinear(np.log10(fith_smooth), euler_params[0], euler_params[1])
    print("Error \u221d h^{0:.3g}".format(euler_params[0]))
    rk4_params, params_covariance = optimize.curve_fit(fitlinear, np.log10(testh), np.log10(rk4errors), p0=[1, 1])
    rk4errorfit = fitlinear(np.log10(fith_smooth), rk4_params[0], rk4_params[1])
    print("Error \u221d h^{0:.3g}".format(rk4_params[0]))

    plt.subplot(1, 3, index)

    plt.scatter(testh, eulererrors, label='Euler')
    plt.plot(fith_smooth, 10 ** eulererrorfit, label='Euler Fit')
    plt.scatter(testh, rk4errors, label='RK4')
    plt.plot(fith_smooth, 10 ** rk4errorfit, label='RK4 Fit')
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('h')
    plt.ylabel('Absolute Error')
    plt.title('n = {0}'.format(j))

plt.show()