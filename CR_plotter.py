"""
    Plotting module.
    Takes care of all plots.
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
import CR_star as star
import CR_config as cf
import CR_laneemden as le

# TODO: [Olivia] Add "useLatex" parameter.
plt.rcParams.update({'font.size': 18, "font.family": "Times New Roman"})


def plotPolytropicDensityRadius(nMin: float = 0.0, nMax: float = 5.5, nStep: float = 0.5,
                                method: str = "rk4", verbose: bool = False) -> None:
    """
    Plots the density-radius plot for a range of polytropic n coefficients.
    Uses ``star.Polytropic(cf.Var.rho).getDensityRadius`` for each n coefficient.
    :param nMin: Minimum range.
    :param nMax: Maximum range. Range will only go up to ``nMax - nStep``
    :param nStep: Step size of the n coefficients.
    :param method: The numerical method that should be used to solve the polytropic solution.
    :param verbose: Print out information on what n is currently worked on.

    **Example:**

    ``plotter.plotPolytropicDensityRadius(nMin = 0.0, nMax = 5.0, nStep = 1.0, method = "rk4")`` will
    plot the density-radius diagram for ``n = 0.0, 1.0, 2.0, 3.0, 4.0`` with Runge-Kutta-4th Order.

    **Thread Safety:**

    NOT SAFE - The function access the global variable ``cf.Var.n``,
    making it impossible to run this in parallel threads.
    """

    # Create the list of n coefficients
    nList = np.arange(nMin, nMax, nStep)

    # TODO: [Marc] Create config setting for figure options
    plt.figure(figsize=(9, 8), dpi=100)

    if verbose:
        print("Plotting density-radius plot for polytropic stars")

    for n in nList:

        # Set the global variable to the different ns
        cf.Var.n = n
        if verbose:
            print("Calculating {0}".format(n))
        r, rho = star.Polytropic(cf.Var.rho).getDensityRadius(method=method, pressure="Non-Relativistic")

        # Divide by the max value to normalize the density and radius
        plt.plot(r / np.max(r), rho / np.max(rho), label="n = {0}".format(n))

    cf.Var.n = 1.5

    plt.xlabel("Dimensionless radius")
    plt.ylabel("Dimensionless density")
    plt.grid()

    # This displays a title for the legend, over all other labels.
    plt.legend(title="n")
    plt.show()


def _cubicInverse(x: float, a: float, b: float, c: float) -> float:
    """
    A fit function for ``plotPolytropicRange``.
    The function is of form `a/((x + b) ** 3) + c`
    :param x: x-coordinate.
    :param a: scaling factor.
    :param b: offset along x-axis.
    :param c: offset along y-axis.
    :return: The y-value of mentioned function.
    """
    return a / ((x + b) ** 3) + c


def printPolytropic(density: float, method: str) -> None:
    """
    Prints out the result of the polytropic differential equation system for just one density.
    Basically, prints out the coordinates of one dot of a mass-radius diagram.
    :param density: Supposed density.
    :param method: The numerical method that should be used to solve the polytropic solution.
    :return: None, but prints out mass radius points for a relativistic
    and non-relativistic pressure function (TOV and Hydro respectively)

    **Example:**
    ``printPolytropic(1e7, "rk4")`` will put out:

        ``For relativistic pressure and rho = 1e7: r = 1e9, m = 5e34``

        ``For relativistic pressure and rho = 1e7: r = 5e9, m = 102e34``

    **Thread safety:**
    NOT TESTED
    """
    star.Polytropic(density) \
        .getMassRadius(method=method, pressure="Relativistic", verbose=True)
    star.Polytropic(density) \
        .getMassRadius(method=method, pressure="Non-Relativistic", verbose=True)


def printWhiteDwarf(density: float, method: str) -> None:
    """
        Prints out the result of the white Dwarf differential equation system for just one density.
        See ``printPolytropic``.
    """

    # Create cubic splines for white dwarfs
    cubicSplineData = np.loadtxt("EoS_Spline_Data.csv", delimiter='\t').transpose()
    rhoEOS = CubicSpline(cubicSplineData[1], cubicSplineData[0], extrapolate=True)
    pEOS = CubicSpline(cubicSplineData[0], cubicSplineData[1], extrapolate=True)

    star.WhiteDwarf(density, pEOS, rhoEOS) \
        .getMassRadius(method=method, pressure="Relativistic", verbose=True)
    star.WhiteDwarf(density, pEOS, rhoEOS) \
        .getMassRadius(method=method, pressure="Non-Relativistic", verbose=True)


def printNeutronStar(density, method):
    """
        Prints out the result of the neutron star differential equation system for just one density.
        See ``printPolytropic``.
    """

    # Create cubic splines for neutron stars
    cubicSplineData = np.loadtxt("wrishikData/Neutron-Star-Structure-master/SLy.txt").transpose()
    rhoEOS = CubicSpline(cubicSplineData[3], cubicSplineData[2], extrapolate=True)
    pEOS = CubicSpline(cubicSplineData[2], cubicSplineData[3], extrapolate=True)

    star.NeutronStar(density, pEOS, rhoEOS) \
        .getMassRadius(method=method, pressure="Relativistic", verbose=True)
    star.NeutronStar(density, pEOS, rhoEOS) \
        .getMassRadius(method=method, pressure="Non-Relativistic", verbose=True)


def plotPolytropicRange(rhoMin: float = 1e7, rhoMax: float = 1e12, rhoNum: int = 30, rhoH: float = 1e5,
                        method: str = "rk4", verbose: bool = False) -> None:
    """
    Plots the mass-radius diagram for a polytropic star for a set range of densities.
    Additionally, it plots and prints out the coefficients for the fit for the Hydro pressure function.
    :param rhoMin: The minimum range of the density.
    :param rhoMax: The maximum range of the density. Range will only go up to ``rhoMax - rhoNum``
    :param rhoNum: The number of density points to be calculated for the graph.
    :param rhoH: The density step size for the differential equation solver.
    :param method: The numerical method that should be used to solve the polytropic solution.
    :param verbose: If additional information should be plotted, like at what calculation the differential equation is currently working on.

    **Thread safety:** NOT TESTED

    :return: Nothing, but makes a mass-radius-plot for the polytropic star.
    """

    rhoValues = np.geomspace(rhoMin, rhoMax, rhoNum)
    relDataValues = np.zeros((rhoValues.size, 2))
    nonRelValues = np.zeros((rhoValues.size, 2))

    for i, rho in enumerate(rhoValues):
        print("Calculation {0}".format(i))
        relDataValues[i, 0], relDataValues[i, 1] = \
            star.Polytropic(rho).getMassRadius(rhoH=rhoH, method=method, verbose=verbose)
        nonRelValues[i, 0], nonRelValues[i, 1] = \
            star.Polytropic(rho).getMassRadius(rhoH=rhoH, method=method, verbose=verbose, pressure="Non-Relativistic")

    # Make the fit for the non-relativistic pressure function
    # TIP: The inline-comment noqa suppresses warnings for that particular line.
    popt, _ = curve_fit(_cubicInverse, nonRelValues[:, 0], nonRelValues[:, 1], p0=[rhoH / 10, 1, 1])  # noqa
    print(popt)

    plt.figure(figsize=(9, 8), dpi=100)
    plt.scatter(nonRelValues[:, 0], nonRelValues[:, 1], color="b", marker="o", label="Hydrostatic")
    plt.scatter(relDataValues[:, 0], relDataValues[:, 1], color="r", marker="v", label="TOV")
    plt.plot(nonRelValues[:, 0], _cubicInverse(nonRelValues[:, 0], *popt))
    plt.xlabel("radius [km]")

    # TODO: [Olivia] Please add the corresponding density units.
    plt.ylabel("Density []")
    plt.grid()
    plt.legend()
    plt.show()


def plotWhiteDwarfRange(rhoMin: float = 1e6, rhoMax: float = 1e14, rhoNum: int = 30, rhoH: float = 1e5,
                        method: str = "rk4", verbose: bool = False) -> None:
    """
        Plots the mass-radius diagram for a white dwarf for a set range of densities.
        Generates ``scipy.CubicSplines`` with the pre-generated data values of the
        dwarfs equation of state described in File ``"EoS_Spline_Data.csv"``.
        :param rhoMin: The minimum range of the density.
        :param rhoMax: The maximum range of the density. Range will only go up to ``rhoMax - rhoNum``
        :param rhoNum: The number of density points to be calculated for the graph.
        :param rhoH: The density step size for the differential equation solver.
        :param method: The numerical method that should be used to solve the white dwarf solution.
        :param verbose: If additional information should be plotted, like at what calculation the differential equation is currently working on.

        **Thread safety:** NOT TESTED

        **Also see:** ``plotPolytropicRange``
        :return: Nothing, but makes a mass-radius-plot for the white dwarf.
    """

    rhoValues = np.geomspace(rhoMin, rhoMax, rhoNum)
    relDataValues = np.zeros((rhoValues.size, 2))
    nonRelValues = np.zeros((rhoValues.size, 2))

    cubicSplineData = np.loadtxt("EoS_Spline_Data.csv", delimiter='\t').transpose()
    rhoEOS = CubicSpline(cubicSplineData[1], cubicSplineData[0], extrapolate=True)
    pEOS = CubicSpline(cubicSplineData[0], cubicSplineData[1], extrapolate=True)

    for i, rho in enumerate(rhoValues):
        print("Calculation {0}".format(i))
        relDataValues[i, 0], relDataValues[i, 1] = \
            star.WhiteDwarf(rho, pEOS, rhoEOS).getMassRadius(rhoH=rhoH, method=method, verbose=verbose)
        nonRelValues[i, 0], nonRelValues[i, 1] = \
            star.WhiteDwarf(rho, pEOS, rhoEOS).getMassRadius(rhoH=rhoH, method=method, verbose=verbose,
                                                             pressure="Non-Relativistic")

    plt.figure(figsize=(9, 8), dpi=100)
    plt.scatter(nonRelValues[:, 0], nonRelValues[:, 1], color="b", marker="o", label="Hydrostatic")
    plt.scatter(relDataValues[:, 0], relDataValues[:, 1], color="r", marker="v", label="TOV")
    whiteDwarfMaxMass = np.max(relDataValues[:, 1])
    plt.axhline(whiteDwarfMaxMass, color="r", linestyle="dashed", label="Chandrasekhar limit")
    plt.xlabel("radius [km]")
    plt.ylabel("mass [$M_{\odot}$]")
    plt.grid()
    plt.legend()
    plt.show()
    print("White Dwarf Max Mass {0}".format(whiteDwarfMaxMass))


def plotNeutronStarRange(rhoMin: float = 2.65e14, rhoMax: float = 1e15, rhoNum: int = 30, rhoH: float = 1e3,
                         method: str = "rk4", verbose: bool = False) -> None:
    """
        Plots the mass-radius diagram for a neutron star for a set range of densities.
        Generates ``scipy.CubicSplines`` with the pre-generated data values of the
        dwarfs equation of state described in File ``"wrishikData/Neutron-Star-Structure-master/SLy.txt"``.
        :param rhoMin: The minimum range of the density.
        :param rhoMax: The maximum range of the density. Range will only go up to ``rhoMax - rhoNum``
        :param rhoNum: The number of density points to be calculated for the graph.
        :param rhoH: The density step size for the differential equation solver.
        :param method: The numerical method that should be used to solve the neutron star solution.
        :param verbose: If additional information should be plotted, like at what calculation the differential equation is currently working on.

        **Thread safety:** NOT TESTED

        **Also see:** ``plotPolytropicRange``
        :return: Nothing, but makes a mass-radius-plot for the white dwarf.
    """

    rhoValues = np.geomspace(rhoMin, 4 * rhoMax, rhoNum)
    relDataValues = np.zeros((rhoValues.size, 2))

    rhoValues2 = np.geomspace(rhoMin, rhoMax, rhoNum)
    nonRelValues = np.zeros((rhoValues2.size, 2))

    cubicSplineData = np.loadtxt("wrishikData/Neutron-Star-Structure-master/SLy.txt").transpose()
    rhoEOS = CubicSpline(cubicSplineData[3], cubicSplineData[2], extrapolate=True)
    pEOS = CubicSpline(cubicSplineData[2], cubicSplineData[3], extrapolate=True)

    for i, rho in enumerate(rhoValues):
        print("Calculation {0}".format(i))
        relDataValues[i, 0], relDataValues[i, 1] = \
            star.NeutronStar(rho, pEOS, rhoEOS).getMassRadius(rhoH=rhoH, method=method, verbose=verbose)

    for i, rho in enumerate(rhoValues2):
        print("Calculation {0}".format(i))
        nonRelValues[i, 0], nonRelValues[i, 1] = \
            star.NeutronStar(rho, pEOS, rhoEOS).getMassRadius(rhoH=4 * rhoH, method=method, verbose=verbose,
                                                              pressure="Non-Relativistic")

    plt.figure(figsize=(9, 8), dpi=100)
    plt.scatter(nonRelValues[:, 0], nonRelValues[:, 1], color="b", marker="o", label="Hydrostatic")
    plt.scatter(relDataValues[:, 0], relDataValues[:, 1], color="r", marker="v", label="TOV")
    plt.xlabel("radius [km]")
    plt.ylabel("mass [$M_{\odot}$]")
    neutronStarMaxMass = np.max(relDataValues[:, 1])
    print("Neutron Star Max Mass {0}".format(neutronStarMaxMass))
    plt.axhline(neutronStarMaxMass, color="r", linestyle="dashed", label="Neutron star mass limit")
    plt.grid()
    plt.legend()
    plt.show()


def plotPolytropicEOS() -> None:
    """
    Prints a density - pressure diagram for the equation os state of polytropic stars.
    :return: Nothing but the plot

    **Thread safety:** NOT NECESSARY
    """

    plt.figure(figsize=(9, 8), dpi=100)
    rho, pressure = star.Polytropic(0).getEOSData()
    plt.plot(rho, pressure, color="k")
    plt.grid()
    plt.xlabel("Density [g/ccm]")

    # TODO: [Olivia] Please add the units for pressure.
    plt.ylabel("Pressure []")
    plt.show()


def plotWhiteDwarfEOS() -> None:
    """
        Prints a density - pressure diagram for the equation os state of white dwarfs.
        :return: Nothing but the plot

        **Thread safety:** NOT NECESSARY

        **Also see:**
        ``plotWhiteDwarfEOS``
    """

    cubicSplineData = np.loadtxt("EoS_Spline_Data.csv", delimiter='\t').transpose()
    rhoEOS = CubicSpline(cubicSplineData[1], cubicSplineData[0], extrapolate=True)
    pEOS = CubicSpline(cubicSplineData[0], cubicSplineData[1], extrapolate=True)

    plt.figure(figsize=(9, 8), dpi=100)
    rho, pressure = star.WhiteDwarf(0, pEOS, rhoEOS).getEOSData()
    plt.plot(rho, pressure, color="k")
    plt.xlabel("Density [g/ccm]")

    # TODO: [Olivia] Please add the units for pressure.
    plt.ylabel("Pressure []")
    plt.show()


def plotNeutronStarEOS():
    """
        Prints a density - pressure diagram for the equation os state of white dwarfs.
        :return: Nothing but the plot

        **Thread safety:** NOT NECESSARY

        **Also see:**
        ``plotWhiteDwarfEOS``
    """
    cubicSplineData = np.loadtxt("wrishikData/Neutron-Star-Structure-master/SLy.txt").transpose()
    rhoEOS = CubicSpline(cubicSplineData[3], cubicSplineData[2], extrapolate=True)
    pEOS = CubicSpline(cubicSplineData[2], cubicSplineData[3], extrapolate=True)

    plt.figure(figsize=(9, 8), dpi=100)
    rho, pressure = star.NeutronStar(0, pEOS, rhoEOS).getEOSData()
    plt.plot(rho, pressure, color="k")

    plt.xlabel("Density [g/ccm]")

    # TODO: [Olivia] Please add the units for pressure.
    plt.ylabel("Pressure []")
    plt.show()


def plotLaneEmden(n: float, analytical: bool, method: str = "rk4") -> None:
    """
    Plots solved solution for a single n.
    Can also plot the analytical solution.
    This function does not plot or show on its own, so if used separately,
    wrap it around with plt.figure() and plt.legend() plt.show()
    :param n: The polytropic coefficient, real positive.
    :param analytical: The analytical solution should be plotted.
    :param method: The numerical method that should be used to solve the lane emden solution.
    :return: Nothing, not even the plot.
    """

    # plt.figure(figsize=(9, 8), dpi=100)
    xi, theta = le.LaneEmden(n).getSolution(stopTime=50, method=method)
    plt.plot(xi, theta, label="${0}$".format(n))
    if analytical:
        plotLaneEmdenAnalytical(n)
    plt.xlabel("Dimensionless radius")
    plt.ylabel("Dimensionless density")
    plt.axhline(0, color='gray', linestyle='dashed')
    plt.grid()
    # plt.legend(title="n")
    # plt.show()


def plotLaneEmdenRange(nMin: float = 0.5, nMax: float = 5.5, nStep: float = 0.5) -> None:
    """
    Plots only the numerical solutions of Lane Emden.
    It does not show the plot, so it has to be followed up with plt.show()
    :param nMin: The lower bound for the polytropic coefficient
    :param nMax: The upper Bound. Range only goes to ``nMax - nStep``
    :param nStep: The step size for n.
    :return: Plot.

    **Warning: This function takes ~10 minutes on good hardware with standard parameters.
    Do not worry about when it takes a bit longer.**
    """

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


def plotLaneEmdenAnalyticalAll() -> None:
    """
    Plots all analytical solutions of Lane-Emden at once.
    This includes n = 0, n = 1, n = 5.
    :return: Nothing but a plot.
    """

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


def plotLaneEmdenAnalytical(n: float):
    """
    Plots only one analytical solution of the Lane Emden Equation.
    The only valid analytical solutions exist for 0, 1, 5.
    If an n is chosen, that is not part of the analytical solution, it will
    print a warning, but will still be part of the plot.
    Additionally, the function has to be wrapped with ``plt.figure(figsize=(9, 8), dpi=100)`` and
    ``plt.legend(title="n"), plt.show()``.

    :param n: The polytropic coefficient.
    :return: Nothing, not even a plot.
    """

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
