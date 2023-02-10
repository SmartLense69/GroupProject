import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
import euler as eu
import rungekutta4lane as rk4l
import rungekutta4 as rk4
import rungekuttafehlberg as rukufe
import secondordermethod as midpoint
import fourthordermethod as fourth
import fifthordermethod as fifth
import rknmethod as rkn

def anader0(xi):
    return -2*xi/6

# gg = np.arange(0, 3, 0.001)
# plt.plot(gg, anader0(gg))
# plt.show()

def anader1(xi):
    return (xi*np.cos(xi) - np.sin(xi))/(xi**2)

def anader5(xi):
    return (-xi/3)*((1+(((xi)**2)/3))**(-3/2))

testh = [1e-1, 5e-1, 1e-2, 5e-2, 1e-3, 5e-3, 1e-4, 5e-4, 5e-5]

# fitting function
def fitlinear(x, a, b):
    return a*x + b

# defining range of values to plot fit function over, sufficient to be smooth line
fith_smooth = np.linspace(min(testh), max(testh), 1000)

plt.figure(figsize=(12, 4), layout='constrained')
plt.suptitle('Error Plots Running to 1 for Derivative')

js = [0, 1, 5]
indices = [1, 2, 3]
analyt = [anader0, anader1, anader5]

for j, index, ana in zip(js, indices, analyt):

    # euler
    eulererrors = np.zeros(len(testh))
    rk4errors = np.zeros(len(testh))
    # rkferrors = np.zeros(len(testh))
    midpointerrors = np.zeros(len(testh))
    fourtherrors = np.zeros(len(testh))
    # fiftherrors = np.zeros(len(testh))
    # rknerrors = np.zeros(len(testh))

    for i in range(0, len(testh)):
        h = testh[i]

        # eulerto1 = eu.euler(3, h, j, eu.phigrad, 0, eu.thetagrad, 1)
        # abserroreulerat1 = np.abs(eulerto1[2][len(eulerto1[2]) - 1] - ana(eulerto1[0][len(eulerto1[0]) - 1]))
        # eulererrors[i] = abserroreulerat1

        # midpointto1 = midpoint.secondorder(3, h, j, midpoint.phigrad, 0, midpoint.thetagrad, 1)
        # abserrormidpointat1 = np.abs(ana(midpointto1[0][len(midpointto1[0]) - 1]) - midpointto1[2][len(midpointto1[2]) - 1])
        # midpointerrors[i] = abserrormidpointat1

        fourthto1 = fourth.fourthorder(3, h, j, fourth.phigrad, 0, fourth.thetagrad, 1)
        abserrorfourthat1 = np.abs(
            ana(fourthto1[0][len(fourthto1[0]) - 1]) - fourthto1[2][len(fourthto1[2]) - 1])
        fourtherrors[i] = abserrorfourthat1

        # rk4to1 = rk4.rungekutta4lane(rk4.phigrad, 3, 0, 1, j, h)
        rk4to1 = rk4.rungekutta4(rk4.phigrad, rk4.thetagrad, rk4.statefunc, 3, 0, 1, j, h)
        abserrorrk4at1 = np.abs(ana(rk4to1[0][-1]) - rk4to1[1][-1])
        rk4errors[i] = abserrorrk4at1

        # fifthto1 = fifth.fifthorder(3, h, j, fifth.phigrad, 0, fifth.thetagrad, 1)
        # abserrorfifthat1 = np.abs(
        #     ana(fifthto1[0][len(fifthto1[0]) - 1]) - fifthto1[1][len(fifthto1[1]) - 1])
        # fiftherrors[i] = abserrorfifthat1

    # print("For n = {0}, the error dependence on h are:".format(j))

    # euler_params, params_covariance = optimize.curve_fit(fitlinear, np.log10(testh), np.log10(eulererrors),
    #                                                      p0=[1, 1])
    # eulererrorfit = fitlinear(np.log10(fith_smooth), euler_params[0], euler_params[1])
    # print("Euler Error \u221d h^{0:.3g}".format(euler_params[0]))

    # midpoint_params, params_covariance = optimize.curve_fit(fitlinear, np.log10(testh), np.log10(midpointerrors),
    #                                                         p0=[1, 1])
    # midpointerrorfit = fitlinear(np.log10(fith_smooth), midpoint_params[0], midpoint_params[1])
    # print("Midpoint Error \u221d h^{0:.3g}".format(midpoint_params[0]))

    fourth_params, params_covariance = optimize.curve_fit(fitlinear, np.log10(testh), np.log10(fourtherrors),
                                                          p0=[1, 1])
    fourtherrorfit = fitlinear(np.log10(fith_smooth), fourth_params[0], fourth_params[1])
    print("3/8-Rule RK4 Error \u221d h^{0:.3g}".format(fourth_params[0]))

    rk4_params, params_covariance = optimize.curve_fit(fitlinear, np.log10(testh), np.log10(rk4errors), p0=[1, 1])
    rk4errorfit = fitlinear(np.log10(fith_smooth), rk4_params[0], rk4_params[1])
    print("RK4 Error \u221d h^{0:.3g}".format(rk4_params[0]))

    print()

    plt.subplot(1, 3, index)

    # plt.scatter(testh, eulererrors, label='Euler', color='purple')

    # plt.plot(fith_smooth, 10 ** eulererrorfit, label='Euler Fit')

    # plt.scatter(testh, midpointerrors, label='Midpoint')
    # plt.plot(fith_smooth, 10 ** midpointerrorfit, label='Midpoint Fit')

    plt.scatter(testh, fourtherrors, label='3/8-Rule RK4')
    plt.plot(fith_smooth, 10 ** fourtherrorfit, label='3/8-Rule Fit')

    plt.scatter(testh, rk4errors, label='Classic RK4')
    plt.plot(fith_smooth, 10 ** rk4errorfit, label='RK4 Fit')

    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('h')
    plt.ylabel('Absolute Error')
    plt.title('n = {0}'.format(j))

plt.show()
