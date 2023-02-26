import sys

from matplotlib import pyplot as plt

import CR_exitcodes as ec
from CR_param import Param as P
import CR_plotter as plotter


class CommandLineInput:

    def _printHelpHint(self):
        print("For more information, type {0} --help or {0} -h.".format(self.argv[0]))

    @staticmethod
    def _printHelp():
        print("There are all valid command line arguments...")
        sys.exit(ec.SUCCESS)

    def _printInvalidArg(self, index):
        print("{0} is not a valid argument in this case.".format(self.argv[index]))
        self._printHelpHint()
        sys.exit(ec.INVALID_ARGUMENT)

    def _checkForArgSyntax(self):
        for i, argv in enumerate(self.argv[1:]):
            if argv[0] != '-':
                self._printInvalidArg(i + 1)

    @staticmethod
    def _printInvalidParam(inputString, param):
        print("{0} is not a valid {1}. A real positive number is required for a {1}.".format(inputString, param))
        sys.exit(ec.INVALID_PARAMETER)

    def _convertToInt(self, inputString, param):
        try:
            intResult = int(inputString)
            if intResult <= 0:
                self._printInvalidParam(inputString=inputString, param=param)
            else:
                return intResult
        except ValueError:
            self._printInvalidParam(inputString=inputString, param=param)

    def _convertToFloat(self, inputString, param):
        try:
            floatResult = float(inputString)
            if floatResult <= 0:
                self._printInvalidParam(inputString=inputString, param=param)
            else:
                return floatResult
        except ValueError:
            self._printInvalidParam(inputString=inputString, param=param)

    def _checkForLaneEmden(self):

        if self.argc > 2:
            if P.N_SINGLE in self.argv[2]:
                self.n = self._convertToFloat(self.argv[2].replace(P.N_SINGLE, ""), param="n-Value")
            if self.argc > 3:
                for i in range(3, self.argc):
                    if P.N_MIN in self.argv[i]:
                        self.nMin = self._convertToFloat(inputString=self.argv[i].replace(P.N_MIN, ""),
                                                         param="n-Value")
                    elif P.N_STEP in self.argv[i]:
                        self.nH = self._convertToFloat(inputString=self.argv[i].replace(P.N_STEP, ""),
                                                       param="n-Value")
                    elif P.N_MAX in self.argv[i]:
                        self.nMax = self._convertToFloat(inputString=self.argv[i].replace(P.N_MAX, ""),
                                                         param="n-Value")
                    elif self.argv[i] == P.RK2[0] or self.argv[i] == P.RK2[1]:
                        self.method = "rk2"
                    elif self.argv[i] == P.RK4[0] or self.argv[i] == P.RK4[1]:
                        self.method = "rk4"
                    elif self.argv[i] == P.RKF[0] or self.argv[i] == P.RKF[1]:
                        self.method = "rkf"
                    elif self.argv[i] == P.EULER[0] or self.argv[i] == P.EULER[1]:
                        self.method = "euler"
                    elif self.argv[i] == P.ANALYTICAL[0] or self.argv[i] == P.ANALYTICAL[1]:
                        self.analytical = True
                    else:
                        self._printInvalidArg(i)
        else:
            print("Missing n-Values for lane-emden.")
            self._printHelpHint()
            sys.exit(ec.MISSING_ARGUMENT)

    def _checkForStar(self):
        if self.argc > 2:
            if P.DENSITY[0] in self.argv[2]:
                self.density = self._convertToFloat(self.argv[2].replace(P.DENSITY[0], ""), param="density")
            elif P.DENSITY[1] in self.argv[2]:
                self.density = self._convertToFloat(self.argv[2].replace(P.DENSITY[1], ""), param="density")
            else:
                for i in range(2, self.argc):
                    if P.DENSITY_MIN[0] in self.argv[i]:
                        self.densityMin = self._convertToFloat(
                            inputString=self.argv[i].replace(P.DENSITY_MIN[0], ""),
                            param="density")
                    elif P.DENSITY_MIN[1] in self.argv[i]:
                        self.densityMin = self._convertToFloat(
                            inputString=self.argv[i].replace(P.DENSITY_MIN[1], ""),
                            param="density")
                    elif P.DENSITY_STEP[0] in self.argv[i]:
                        self.densityStep = self._convertToFloat(
                            inputString=self.argv[i].replace(P.DENSITY_STEP[0], ""),
                            param="density")
                    elif P.DENSITY_STEP[1] in self.argv[i]:
                        self.densityStep = self._convertToFloat(
                            inputString=self.argv[i].replace(P.DENSITY_STEP[1], ""),
                            param="density")
                    elif P.DENSITY_MAX[0] in self.argv[i]:
                        self.densityMax = self._convertToFloat(
                            inputString=self.argv[i].replace(P.DENSITY_MAX[0], ""),
                            param="density")
                    elif P.DENSITY_MAX[1] in self.argv[i]:
                        self.densityMax = self._convertToFloat(
                            inputString=self.argv[i].replace(P.DENSITY_MAX[1], ""),
                            param="density")
                    elif P.DENSITY_NUM[0] in self.argv[i]:
                        self.densityNum = self._convertToInt(
                            inputString=self.argv[i].replace(P.DENSITY_NUM[0], ""),
                            param="density")
                    elif P.DENSITY_NUM[1] in self.argv[i]:
                        self.densityNum = self._convertToInt(
                            inputString=self.argv[i].replace(P.DENSITY_NUM[1], ""),
                            param="density")
                    elif P.N_MIN in self.argv[i]:
                        self.nMin = self._convertToFloat(inputString=self.argv[i].replace(P.N_MIN, ""),
                                                         param="n-Value")
                    elif P.N_STEP in self.argv[i]:
                        self.nH = self._convertToFloat(inputString=self.argv[i].replace(P.N_STEP, ""),
                                                       param="n-Value")
                    elif P.N_MAX in self.argv[i]:
                        self.nMax = self._convertToFloat(inputString=self.argv[i].replace(P.N_MAX, ""),
                                                         param="n-Value")
                    elif self.argv[i] == P.RK2[0] or self.argv[i] == P.RK2[1]:
                        self.method = "rk2"
                    elif self.argv[i] == P.RK4[0] or self.argv[i] == P.RK4[1]:
                        self.method = "rk4"
                    elif self.argv[i] == P.RKF[0] or self.argv[i] == P.RKF[1]:
                        self.method = "rkf"
                    elif self.argv[i] == P.EULER[0] or self.argv[i] == P.EULER[1]:
                        self.method = "euler"
                    elif self.argv[i] == P.SHOW_EOS[0] or self.argv[i] == P.SHOW_EOS[1]:
                        self.showEOS = True
                    else:
                        self._printInvalidArg(i)
        else:
            print("Missing densities for given star.")
            self._printHelpHint()
            sys.exit(ec.MISSING_ARGUMENT)

    def _sanitizeInput(self):
        nProvided = False
        if self.method is None:
            print("Please provide a numerical method.")
            sys.exit(ec.MISSING_ARGUMENT)
        else:
            if self.method not in ["rk2", "rk4", "rkf", "euler"]:
                print("The numerical method {0} is not supported.".format(self.method))
        if self.laneEmden is not None:
            if self.starType is not None:
                print("Please specify either what star type or lane emden solution is needed.")
                sys.exit(ec.INVALID_ARGUMENT)
            if self.n is not None:
                if self.nH is not None or self.nMin is not None or self.nMax is not None:
                    print("Either specify a single n solution, or give a range.")
                    sys.exit(ec.INVALID_PARAMETER)
            elif self.nH is None or self.nMin is None or self.nMax is None:
                print("Please provide a range for n of the form: min, step, max")
                sys.exit(ec.MISSING_ARGUMENT)
        elif self.starType is not None:
            if self.laneEmden is not None:
                print("Please specify either what star type or lane emden solution is needed.")
                sys.exit(ec.INVALID_ARGUMENT)
            if self.starType not in ["polytropic", "whiteDwarf", "neutronStar"]:
                print("The star type {0} is not supported".format(self.starType))
            if self.starType == "polytropic":
                if self.nMin is not None or self.nMax is not None or self.nH is not None:
                    # I am not negating this statement.
                    if self.nMin is not None and self.nMax is not None and self.nH is not None:
                        nProvided = True
                    else:
                        print("Please provide a range for n of the form: min, step, max")
                        sys.exit(ec.MISSING_ARGUMENT)
            if self.density is not None:
                if self.densityNum is not None or self.densityMin \
                        is not None or self.densityMax is not None or self.densityNum:
                    print("Either specify a single density, or give a range.")
                    sys.exit(ec.INVALID_PARAMETER)
            elif (self.densityNum is None or self.densityMin
                    is None or self.densityMax is None or self.densityNum is None) and nProvided is False:
                print("Please provide a range for density of the form: min, step, max, num")
                sys.exit(ec.MISSING_ARGUMENT)

        else:
            print("Please specify first what star type or lane emden solution is needed.")
            sys.exit(ec.MISSING_ARGUMENT)

    def execute(self):
        if self.laneEmden is not None:
            plt.figure(figsize=(9, 8), dpi=100)
            if self.n is not None:
                if self.analytical is not None:
                    plotter.plotLaneEmden(n=self.n, analytical=self.analytical, method=self.method)
                else:
                    plotter.plotLaneEmden(n=self.n, analytical=False, method=self.method)
            else:
                if self.analytical is not None:
                    plotter.plotLaneEmdenAnalyticalAll()
                else:
                    plotter.plotLaneEmdenRange(self.nMin, self.nMax, self.nH)
            plt.grid()
            plt.legend(title="$n$")
            plt.show()
        elif self.starType is not None:
            if self.showEOS is not None:
                match self.starType:
                    case "polytropic":
                        plotter.plotPolytropicEOS()
                    case "whiteDwarf":
                        plotter.plotWhiteDwarfEOS()
                    case "neutronStar":
                        plotter.plotNeutronStarEOS()
            if self.density is not None:
                match self.starType:
                    case "polytropic":
                        plotter.printPolytropic(self.density, method=self.method)
                    case "whiteDwarf":
                        plotter.printPolytropic(self.density, method=self.method)
                    case "neutronStar":
                        plotter.printPolytropic(self.density, method=self.method)
            if self.starType == "polytropic" and \
                (self.nMin is not None and self.nMax is not None and self.nH is not None):
                plotter.plotPolytropicDensityRadius(nMin=self.nMin, nMax=self.nMax, nStep=self.nH, method=self.method)
            else:
                match self.starType:
                    case "polytropic":
                        plotter.plotPolytropicRange(rhoMin=self.densityMin, rhoMax=self.densityMax,
                                                    rhoNum=self.densityNum, rhoH=self.densityStep, method=self.method)
                    case "whiteDwarf":
                        plotter.plotWhiteDwarfRange(self.densityMin, self.densityMax, self.densityNum, method=self.method)
                    case "neutronStar":
                        plotter.plotNeutronStarRange(self.densityMin, self.densityMax, self.densityNum, method=self.method)

    def __init__(self, argcP: int, argvP: list[str]):

        self.argc: int = argcP
        self.argv: list[str] = argvP
        self.n = None
        self.nMin = None
        self.nMax = None
        self.nH = None
        self.density = None
        self.densityMin = None
        self.densityMax = None
        self.densityStep = None
        self.densityNum = None
        self.method = None
        self.starType = None
        self.laneEmden = None
        self.analytical = None
        self.showEOS = None

        if self.argc == 1:
            print("No arguments were given")
            sys.exit(ec.NO_ARGUMENTS)
        else:
            self._checkForArgSyntax()
            if self.argv[1] == P.LANE_EMDEN[0] or self.argv[1] == P.LANE_EMDEN[1]:
                self.laneEmden = True
                self._checkForLaneEmden()
            elif self.argv[1] == P.POLYTROPIC[0] or self.argv[1] == P.POLYTROPIC[1]:
                self.starType = "polytropic"
                self._checkForStar()
            elif self.argv[1] == P.WHITE_DWARF[0] or self.argv[1] == P.WHITE_DWARF[1]:
                self.starType = "whiteDwarf"
                self._checkForStar()
            elif self.argv[1] == P.NEUTRON_STAR[0] or self.argv[1] == P.NEUTRON_STAR[1]:
                self.starType = "neutronStar"
                self._checkForStar()
            elif self.argv[1] == P.HELP[0] or self.argv[1] == P.HELP[1]:
                self._printHelp()
            else:
                self._printInvalidArg(1)
            self._sanitizeInput()

    def __str__(self):
        output = "Input Parameters\n" \
                 "================\n\n" \
                 "Lane-Emden:\t{0}\n" \
                 "N - single:\t{1}\n" \
                 "N - min:\t{2}\n" \
                 "N - step:\t{3}\n" \
                 "N - max:\t{4}\n\n" \
                 "Star-Type:\t{5}\n" \
                 "rho:\t\t{6}\n" \
                 "rho - min:\t{7}\n" \
                 "rho - step:\t{8}\n" \
                 "rho - max:\t{9}\n" \
                 "Method: \t{10}\n".format(self.laneEmden, self.n, self.nMin, self.nH, self.nMax,
                                           self.starType, self.density, self.densityMin, self.densityNum,
                                           self.densityMax,
                                           self.method)
        return output
