import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import CubicSpline
import CR_star as star


def plotPolytropic():

    rhoValues = np.geomspace(1e7, 1e12, 30)
    relDataValues = np.zeros((rhoValues.size, 2))
    nonRelValues = np.zeros((rhoValues.size, 2))

    for i, rho in enumerate(rhoValues):
        print("Calculation {0}".format(i))
        relDataValues[i, 0], relDataValues[i, 1] = star.Polytropic(rho).getMassRadius()
        nonRelValues[i, 0], nonRelValues[i, 1] = star.Polytropic(rho).getMassRadius(
            pressure="Non-Relativistic")

    plt.rcParams.update({'font.size': 18, "font.family": "Times New Roman"})
    plt.figure(figsize=(9, 8), dpi=100)
    plt.scatter(relDataValues[:, 0], relDataValues[:, 1], color="r", marker="v", label="TOV")
    plt.scatter(nonRelValues[:, 0], nonRelValues[:, 1], color="b", marker="o", label="Hydrostatic")
    plt.xlabel("radius [km]")
    plt.ylabel("mass [$M_{\odot}$]")
    plt.grid()
    plt.legend()
    plt.show()


def plotWhiteDwarf():

    rhoValues = np.geomspace(1e6, 1e14, 30)
    relDataValues = np.zeros((rhoValues.size, 2))
    nonRelValues = np.zeros((rhoValues.size, 2))

    cubicSplineData = np.loadtxt("EoS_Spline_Data.csv", delimiter='\t').transpose()
    rhoEOS = CubicSpline(cubicSplineData[1], cubicSplineData[0], extrapolate=True)
    pEOS = CubicSpline(cubicSplineData[0], cubicSplineData[1], extrapolate=True)

    for i, rho in enumerate(rhoValues):
        print("Calculation {0}".format(i))
        relDataValues[i, 0], relDataValues[i, 1] = star.WhiteDwarf(rho, pEOS, rhoEOS).getMassRadius()
        nonRelValues[i, 0], nonRelValues[i, 1] = star.WhiteDwarf(rho, pEOS, rhoEOS).getMassRadius(pressure="Non-Relativistic")

    plt.rcParams.update({'font.size': 18, "font.family": "Times New Roman"})
    plt.figure(figsize=(9, 8), dpi=100)
    plt.scatter(relDataValues[:, 0], relDataValues[:, 1], color="r", marker="v", label="TOV")
    plt.scatter(nonRelValues[:, 0], nonRelValues[:, 1], color="b", marker="o", label="Hydrostatic")
    plt.xlabel("radius [km]")
    plt.ylabel("mass [$M_{\odot}$]")
    plt.grid()
    plt.legend()
    plt.show()


def plotNeutronStar():

    rhoValues = np.geomspace(2.65e14, 4e15, 30)
    relDataValues = np.zeros((rhoValues.size, 2))

    rhoValues2 = np.geomspace(2.65e14, 1e15, 30)
    nonRelValues = np.zeros((rhoValues2.size, 2))

    cubicSplineData = np.loadtxt("wrishikData/Neutron-Star-Structure-master/SLy.txt").transpose()
    rhoEOS = CubicSpline(cubicSplineData[3], cubicSplineData[2], extrapolate=True)
    pEOS = CubicSpline(cubicSplineData[2], cubicSplineData[3], extrapolate=True)

    for i, rho in enumerate(rhoValues):
        print("Calculation {0}".format(i))
        relDataValues[i, 0], relDataValues[i, 1] = star.WhiteDwarf(rho, pEOS, rhoEOS).getMassRadius(rhoH=1e3)

    for i, rho in enumerate(rhoValues2):
        print("Calculation {0}".format(i))
        nonRelValues[i, 0], nonRelValues[i, 1] =\
            star.WhiteDwarf(rho, pEOS, rhoEOS).getMassRadius(rhoH=4e3, pressure="Non-Relativistic")

    plt.rcParams.update({'font.size': 18, "font.family": "Times New Roman"})
    plt.figure(figsize=(9, 8), dpi=100)
    plt.scatter(relDataValues[:, 0], relDataValues[:, 1], color="r", marker="v", label="TOV")
    plt.scatter(nonRelValues[:, 0], nonRelValues[:, 1], color="b", marker="o", label="Hydrostatic")
    plt.xlabel("radius [km]")
    plt.ylabel("mass [$M_{\odot}$]")
    plt.grid()
    plt.legend()
    plt.show()