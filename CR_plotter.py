import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
import CR_star as star
import CR_config as cf
import CR_laneemden as le

plt.rcParams.update({'font.size': 18, "font.family": "Times New Roman"})


def plotPolytropicRhoRadius(nMin=0, nMax=5.5, nStep=0.5):

    nList = np.arange(nMin, nMax, nStep)

    plt.rcParams.update({'font.size': 18, "font.family": "Times New Roman"})
    plt.figure(figsize=(9, 8), dpi=100)

    for n in nList:
        cf.Var.n = n
        r, rho = star.Polytropic(cf.Var.rho).getDensityRadius(method="rk4", pressure="Non-Relativistic")
        plt.plot(r/np.max(r), rho/np.max(rho), label="n = {0}".format(n))

    cf.Var.n = 1.5

    plt.xlabel("radius [km]")
    plt.ylabel("Dimensionless density")
    plt.grid()
    plt.legend()
    plt.show()


def _cubicInverse(x, a, b, c):
    return a/((x + b) ** 3) + c


def printPolytropic(density, method):
    star.Polytropic(density) \
        .getMassRadius(method=method, pressure="Relativistic", verbose=True)
    star.Polytropic(density) \
        .getMassRadius(method=method, pressure="Non-Relativistic", verbose=True)


def printWhiteDwarf(density, method):

    cubicSplineData = np.loadtxt("EoS_Spline_Data.csv", delimiter='\t').transpose()
    rhoEOS = CubicSpline(cubicSplineData[1], cubicSplineData[0], extrapolate=True)
    pEOS = CubicSpline(cubicSplineData[0], cubicSplineData[1], extrapolate=True)

    star.WhiteDwarf(density, pEOS, rhoEOS) \
        .getMassRadius(method=method, pressure="Relativistic", verbose=True)
    star.WhiteDwarf(density, pEOS, rhoEOS) \
        .getMassRadius(method=method, pressure="Non-Relativistic", verbose=True)


def printNeutronStar(density, method):
    cubicSplineData = np.loadtxt("wrishikData/Neutron-Star-Structure-master/SLy.txt").transpose()
    rhoEOS = CubicSpline(cubicSplineData[3], cubicSplineData[2], extrapolate=True)
    pEOS = CubicSpline(cubicSplineData[2], cubicSplineData[3], extrapolate=True)

    star.NeutronStar(density, pEOS, rhoEOS) \
        .getMassRadius(method=method, pressure="Relativistic", verbose=True)
    star.NeutronStar(density, pEOS, rhoEOS) \
        .getMassRadius(method=method, pressure="Non-Relativistic", verbose=True)


def plotPolytropicRange(rhoMin=1e7, rhoMax=1e12, rhoNum=30, rhoH=1e5, method="rk4", verbose=False):

    rhoValues = np.geomspace(rhoMin, rhoMax, rhoNum)

    relDataValues = np.zeros((rhoValues.size, 2))
    nonRelValues = np.zeros((rhoValues.size, 2))

    for i, rho in enumerate(rhoValues):
        print("Calculation {0}".format(i))
        relDataValues[i, 0], relDataValues[i, 1] =\
            star.Polytropic(rho).getMassRadius(rhoH=rhoH, method=method, verbose=verbose)
        nonRelValues[i, 0], nonRelValues[i, 1] =\
            star.Polytropic(rho).getMassRadius(rhoH=rhoH, method=method, verbose=verbose, pressure="Non-Relativistic")

    popt, _ = curve_fit(_cubicInverse, nonRelValues[:, 0], nonRelValues[:, 1], p0=[rhoH/10, 1, 1])
    print(popt)

    plt.figure(figsize=(9, 8), dpi=100)
    plt.scatter(nonRelValues[:, 0], nonRelValues[:, 1], color="b", marker="o", label="Hydrostatic")
    plt.scatter(relDataValues[:, 0], relDataValues[:, 1], color="r", marker="v", label="TOV")
    plt.plot(nonRelValues[:, 0], _cubicInverse(nonRelValues[:, 0], *popt))
    plt.xlabel("radius [km]")
    plt.ylabel("Density [Pascal]")
    plt.grid()
    plt.legend()
    plt.show()


def plotWhiteDwarf(rhoMin=1e6, rhoMax=1e14, rhoNum=30, rhoH=1e5, method="rk4", verbose=False):

    rhoValues = np.geomspace(rhoMin, rhoMax, rhoNum)
    relDataValues = np.zeros((rhoValues.size, 2))
    nonRelValues = np.zeros((rhoValues.size, 2))

    cubicSplineData = np.loadtxt("EoS_Spline_Data.csv", delimiter='\t').transpose()
    rhoEOS = CubicSpline(cubicSplineData[1], cubicSplineData[0], extrapolate=True)
    pEOS = CubicSpline(cubicSplineData[0], cubicSplineData[1], extrapolate=True)

    for i, rho in enumerate(rhoValues):
        print("Calculation {0}".format(i))
        relDataValues[i, 0], relDataValues[i, 1] =\
            star.WhiteDwarf(rho, pEOS, rhoEOS).getMassRadius(rhoH=rhoH, method=method, verbose=verbose)
        nonRelValues[i, 0], nonRelValues[i, 1] =\
            star.WhiteDwarf(rho, pEOS, rhoEOS).getMassRadius(rhoH=rhoH, method=method, verbose=verbose,
                                                             pressure="Non-Relativistic")

    plt.figure(figsize=(9, 8), dpi=100)
    plt.scatter(nonRelValues[:, 0], nonRelValues[:, 1], color="b", marker="o", label="Hydrostatic")
    plt.scatter(relDataValues[:, 0], relDataValues[:, 1], color="r", marker="v", label="TOV")
    plt.xlabel("radius [km]")
    plt.ylabel("mass [$M_{\odot}$]")
    plt.grid()
    plt.legend()
    plt.show()


def plotNeutronStar(rhoMin=2.65e14, rhoMax=1e15, rhoNum=30, rhoH=1e3,  method="rk4", verbose=False):

    rhoValues = np.geomspace(rhoMin, 4 * rhoMax, rhoNum)
    relDataValues = np.zeros((rhoValues.size, 2))

    rhoValues2 = np.geomspace(rhoMin, rhoMax, rhoNum)
    nonRelValues = np.zeros((rhoValues2.size, 2))

    cubicSplineData = np.loadtxt("wrishikData/Neutron-Star-Structure-master/SLy.txt").transpose()
    rhoEOS = CubicSpline(cubicSplineData[3], cubicSplineData[2], extrapolate=True)
    pEOS = CubicSpline(cubicSplineData[2], cubicSplineData[3], extrapolate=True)

    for i, rho in enumerate(rhoValues):
        print("Calculation {0}".format(i))
        relDataValues[i, 0], relDataValues[i, 1] =\
            star.NeutronStar(rho, pEOS, rhoEOS).getMassRadius(rhoH=rhoH, method=method, verbose=verbose)

    for i, rho in enumerate(rhoValues2):
        print("Calculation {0}".format(i))
        nonRelValues[i, 0], nonRelValues[i, 1] =\
            star.NeutronStar(rho, pEOS, rhoEOS).getMassRadius(rhoH=4 * rhoH, method=method, verbose=verbose,
                                                              pressure="Non-Relativistic")

    plt.figure(figsize=(9, 8), dpi=100)
    plt.scatter(nonRelValues[:, 0], nonRelValues[:, 1], color="b", marker="o", label="Hydrostatic")
    plt.scatter(relDataValues[:, 0], relDataValues[:, 1], color="r", marker="v", label="TOV")
    plt.xlabel("radius [km]")
    plt.ylabel("mass [$M_{\odot}$]")
    plt.grid()
    plt.legend()
    plt.show()


def plotPolytropicEOS():

    plt.figure(figsize=(9, 8), dpi=100)
    rho, pressure = star.Polytropic(0).getEOSData()
    plt.plot(rho, pressure, color="k")
    plt.grid()
    plt.xlabel("Density [g/ccm]")
    plt.ylabel("Pressure")
    plt.show()


def plotWhiteDwarfEOS():

    cubicSplineData = np.loadtxt("EoS_Spline_Data.csv", delimiter='\t').transpose()
    rhoEOS = CubicSpline(cubicSplineData[1], cubicSplineData[0], extrapolate=True)
    pEOS = CubicSpline(cubicSplineData[0], cubicSplineData[1], extrapolate=True)

    plt.figure(figsize=(9, 8), dpi=100)
    rho, pressure = star.WhiteDwarf(0, pEOS, rhoEOS).getEOSData()
    plt.plot(rho, pressure, color="k")
    plt.xlabel("Density [g/ccm]")
    plt.ylabel("Pressure")
    plt.show()


def plotNeutronStarEOS():

    cubicSplineData = np.loadtxt("wrishikData/Neutron-Star-Structure-master/SLy.txt").transpose()
    rhoEOS = CubicSpline(cubicSplineData[3], cubicSplineData[2], extrapolate=True)
    pEOS = CubicSpline(cubicSplineData[2], cubicSplineData[3], extrapolate=True)

    plt.figure(figsize=(9, 8), dpi=100)
    rho, pressure = star.NeutronStar(0, pEOS, rhoEOS).getEOSData()
    plt.plot(rho, pressure, color="k")
    plt.xlabel("Density [g/ccm]")
    plt.ylabel("Pressure")
    plt.show()


def plotLaneEmden(n, analytical):

    # plt.figure(figsize=(9, 8), dpi=100)
    xi, theta = le.LaneEmden(n).getSolution(stopTime=50)
    plt.plot(xi, theta, label="${0}$".format(n))
    if analytical:
        plotLaneEmdenAnalytical(n)
    plt.xlabel("Dimensionless radius")
    plt.ylabel("Dimensionless density")
    plt.axhline(0, color='gray', linestyle='dashed')
    plt.grid()
    # plt.legend(title="n")
    # plt.show()


def plotLaneEmdenRange(nMin=0.5, nMax=5.5, nStep=0.5):

    plt.figure(figsize=(9, 8), dpi=100)

    for n in np.arange(nMin, nMax, nStep):
        xi, theta = le.LaneEmden(n).getSolution(stopTime=50)
        plt.plot(xi, theta, label="${0}$".format(n))
    plt.xlabel("Dimensionless radius")
    plt.ylabel("Dimensionless density")
    plt.axhline(0, color='gray', linestyle='dashed')
    plt.grid()
    plt.legend(title="n")
    # plt.show()


def plotLaneEmdenAnalyticalAll():

    plt.figure(figsize=(9, 8), dpi=100)
    for n in [0, 1, 5]:
        xi, theta = le.LaneEmden(n).getAnalyticalSolution()
        plt.plot(xi, theta, label="{0}".format(n))
    plt.xlabel("Dimensionless radius")
    plt.ylabel("Dimensionless density")
    plt.ylim(-0.25, 1.05)
    plt.axhline(0, color='gray', linestyle='dashed')
    plt.grid()
    plt.legend(title="n")
    # plt.show()


def plotLaneEmdenAnalytical(n):
    # plt.figure(figsize=(9, 8), dpi=100)
    match n:
        case 0.0:
            xi, theta = le.LaneEmden(n).getAnalyticalSolution()
            plt.plot(xi, theta, label="${0}$ Analy.".format(n))
        case 1.0:
            xi, theta = le.LaneEmden(n).getAnalyticalSolution()
            plt.plot(xi, theta, label="${0}$ Analy.".format(n))
        case 5.0:
            xi, theta = le.LaneEmden(n).getAnalyticalSolution()
            plt.plot(xi, theta, label="${0}$ Analy.".format(n))
        case other:
            print("Warning: There are no analytical solution for {0}".format(n))
    plt.xlabel("Dimensionless radius")
    plt.ylabel("Dimensionless density")
    plt.ylim(-0.25, 1.05)
    plt.axhline(0, color='gray', linestyle='dashed')
    plt.grid()
    # plt.legend(title="n")
    # plt.show()

