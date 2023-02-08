import warnings as wr
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import CubicSpline

import config_CR as cf

# np.seterr(all='raise')

_RHO = 0
_P = 0


class InputVariableDuplicateError(Exception):
    """
    Thrown when a differential equation has input variables
    that share the same parameter name
    """

    def __init__(self, nonUniqueVariableName: str):
        self._nonUniqueVariableName = nonUniqueVariableName
        self._message = "A differential equation should have no input variable names that share the same name." \
                        "The input variable name in question is: "
        super().__init__(self._message + self._nonUniqueVariableName)


class InputVariableNumberMismatchError(Exception):
    """
    Thrown when a differential equation was called with
    too many or to little input arguments.
    """

    def __init__(self, sizeDictionary: int, sizeInputVariables: int):
        self._lenInputDictionary = sizeDictionary
        self._lenInputVariables = sizeInputVariables
        self._message = "The differential equation requires " \
                        "{0} inputs but {1} were given"\
            .format(self._lenInputVariables, self._lenInputDictionary)
        super().__init__(self._message)


class SortError(Exception):
    """
    Thrown when the list of parameters are not sorted.
    """

    def __init__(self, flag: int):
        match flag:
            case 0:
                super().__init__("The input variable names list is not sorted")
            case 1:
                super().__init__("The differential equations are not correctly "
                                 "sorted by their calculation order")


class DifferentialEquation:
    """
    A differential equation is a mathematical function,
    where the result is a gradient of an output variable is dependent
    on a series of input variables, including a variable that resembles time.
    """

    def __checkInputDuplicates(self) -> None:
        """
        Checks if there are duplicates in the input variables,
        when initializing the class.
        :raises InputVariableDuplicateError: when there are duplicates in input variables
        """
        for index in range(0, self.inputVarSize - 1):
            currentInputVarName: str = self.inputVarNames[index]
            for compareIndex in range(index+1, self.inputVarSize):
                if currentInputVarName == self.inputVarNames[compareIndex]:
                    raise InputVariableDuplicateError(currentInputVarName)

    def __checkInputSorting(self) -> None:
        """
        Checks if the input variable names are sorted alphabetically
        They have to be alphabetically sorted, as of 04/02/2023,
        this class assigns input values via a dictionary, that is sorted alphabetically.
        TODO: Implement method so input variables do not have to be sorted alphabetically.
        TODO: The sort() method sorts after ASCII occurrence, not after the actual alphabet. USE ONLY SMALL LETTERS
        """
        sortedInput = self.inputVarNames.copy()
        sortedInput.sort()
        for sortedInputItem, inputVarNamesItem in zip(sortedInput, self.inputVarNames):
            if sortedInputItem is not inputVarNamesItem:
                raise SortError(0)

    def call(self, inputValues: dict) -> float:
        """
        Calls the differential equation with values put in.
        :param inputValues: Dictionary of the input values. Syntax is: {a: value1, b: value2}
        :return: Output of the differential equation.
        """
        sortedInputValues = dict(sorted(inputValues.items()))
        sizeDictionary = len(sortedInputValues)

        # Check if the input variables has the same size of the dictionary keyword list
        if self.inputVarSize != sizeDictionary:
            raise InputVariableNumberMismatchError(sizeDictionary, self.inputVarSize)

        # Whenever a numpy array is created, always use np.float128
        numpyInputValues = np.zeros(sizeDictionary, dtype=np.float128)
        for index in range(0, sizeDictionary):
            numpyInputValues[index] = list(sortedInputValues.values())[index]
        return self.functionCall(numpyInputValues)

    def __init__(self, inputVariables: list[str], outputVariable: str, function: callable(np.ndarray),
                 initialCondition: float, inputVarTimeIndex: int):
        """
        Creates a differential equation object.

        Syntax is as follows:
        diff = DifferentialEquation(["input1", "input2"], "output", function, initial, index)

        :param inputVariables: Input variables as a list of strings. Syntax: ["a","b","c"].
        :param outputVariable: For what is this differential equation solving the gradient for.
        :param function: The function in which takes a numpy array of inputs and gives an output.
        :param initialCondition: The initial condition of the differential equation at t = 0.
        :param inputVarTimeIndex: The index of the time dependent variable in inputVariables.
        """
        self.inputVarNames: list[str] = inputVariables
        self.inputVarSize: int = len(self.inputVarNames)
        self.outputVarName: str = outputVariable
        self.functionCall: callable(np.ndarray) = function
        self.initialCondition: float = initialCondition
        self.timeVar: str = inputVariables[inputVarTimeIndex]
        self.__checkInputDuplicates()
        self.__checkInputSorting()


class DifferentialEquationSystem:
    """
    A set of differential equations.
    """

    def __checkSorting(self) -> None:
        """
        **@DEPRECATED: Function is probably not necessary,
        but results after deletion was not tested.**

        Checks if the give differential equations would be computable,
        given in the order they have been passed to the class constructor.
        :raises SortError:
            When a differential equation, which input is computed/output by another differential
            equation after the current differential equation. Example: a = func(b, c), followed by b = func(d)
            -> Throws SortError.
            Should have been: b = func(d) and then a = func(b,c)
        """

        # This function is deprecated
        wr.warn("Function __checkSorting(self) is probably not necessary, but results after deletion was not tested.",
                DeprecationWarning)
        for diffIndex, diff in enumerate(self.listDiffs):
            for compareIndex, compareDiff in enumerate(self.listDiffs):

                # This if statement avoids comparing a differential equation to itself
                if diff is not compareDiff:
                    for inputVars in diff.inputVarNames:

                        # If our comparison differential equation has an output
                        # which is part of the inputs of the current differential equation, and
                        # it comes after this current differential equation, there is no way
                        # to compute the comparison differential equation, since the inputs of the
                        # current differential equation will never be met.
                        if compareDiff.outputVarName is inputVars and compareIndex > diffIndex:
                            raise SortError(1)

    def __unifyVariables(self) -> None:
        """
        Gathers all variables mentioned in the differential equations,
        including input, output and time dependent variables,
        sorts out any duplicates, and saves them in seperate lists.
        """

        # Unify Input Variables
        self.listInput: list[str] = []
        for diff in self.listDiffs:
            for inputVar in diff.inputVarNames:
                self.listInput.append(inputVar)

        # Hack: A dict can only have unique variables.
        # Create one to get rid of duplicates,
        # and make a list out of the dict
        self.listInput = list(dict.fromkeys(self.listInput))

        # Unify time dependent variables
        self.listTime: list[str] = []
        for diff in self.listDiffs:
            self.listTime.append(diff.timeVar)
        self.listTime = list(dict.fromkeys(self.listTime))

        # Unify output variables
        self.listOutput: list[str] = []
        for diff in self.listDiffs:
            self.listOutput.append(diff.outputVarName)
        self.listOutput = list(dict.fromkeys(self.listOutput))

        self.listAllVariables: list[str] = []
        for inputVar in self.listInput:
            self.listAllVariables.append(inputVar)
        for outputVar in self.listOutput:
            self.listAllVariables.append(outputVar)
        self.listAllVariables = list(dict.fromkeys(self.listAllVariables))

    def __init__(self, listDifferentials: list[DifferentialEquation]):
        """
        Create a differential equation system.
        :param listDifferentials: List of differential equation objects
        """
        self.listDiffs: list[DifferentialEquation] = listDifferentials
        # self.__checkSorting()
        self.__unifyVariables()


class DifferentialSolver:
    """
    A class that has several methods to solve a differential equation system.
    """

    varDict: dict[str, np.ndarray]
    """
    A dictionary to keep track of all involved variables.
    
    {key: variable name | value: data (numpy array)}
    """

    thresholdDict: dict[str, float]
    """
    Dictionary for all threshold conditions
    
    {key: variable name | value: threshold value}
    """

    def __resetVarDict(self) -> None:
        """
        Resets the variable dictionary to an ever-increasing numpy array.
        """
        for variable in self.diffSystem.listAllVariables:
            self.varDict[variable] = self.steps.copy()
        for outputVariable, diff in zip(self.diffSystem.listOutput, self.diffSystem.listDiffs):
            self.varDict[outputVariable][0] = diff.initialCondition

    def __init__(self, diffSystem: DifferentialEquationSystem, stepSize: float, stopTime: float):
        """
        Creates a differential equation solver object.
        :param diffSystem: A set of differential equations
        :param stepSize: Representative of h
        :param stopTime: How far should the solver go.
        """

        self.varDict = {}
        self.thresholdDict = {}
        self.diffSystem = diffSystem
        self.steps = np.arange(1, stopTime+stepSize, stepSize)
        self.stepNum = self.steps.size
        """
        The number of elements in the value array in varDict. 
        """

        self.stepSize = stepSize
        """
        The size of one step, representive of variable *h*
        """

        self.__resetVarDict()

    def euler(self) -> None:

        # Reset the variable dictionary, so we don't use a solution
        # computed from a different iteration/method
        self.__resetVarDict()

        # Iterate through the value arrays
        for i in range(1, self.stepNum):

            # Iterate through the differential equations
            for diff in self.diffSystem.listDiffs:

                # Because all functions take input variables as a numpy array,
                # we append the results in a list, because numpy has no append.
                inputList = []
                for inputVariable in diff.inputVarNames:
                    inputList.append(self.varDict.get(inputVariable)[i - 1])

                # Convert to an numpy compatible list
                inputList = np.asarray(inputList)

                # Get the output variable of the current differential equation
                outputArray = self.varDict.get(diff.outputVarName)

                # That's Euler for you
                outputArray[i] = outputArray[i - 1] + self.stepSize * diff.functionCall(inputList)

                if self.thresholdDict.get(diff.outputVarName) is not None:
                    if outputArray[i] < self.thresholdDict.get(diff.outputVarName):
                        self.varDict.update({diff.outputVarName: outputArray})
                        for variable in self.varDict:
                            self.varDict.update({variable: self.varDict.get(variable)[:i-1]})
                        return

                # Update the entry for the output variable
                self.varDict.update({diff.outputVarName: outputArray})

    def addThreshold(self, thresholdDictionary: dict[str, float]):
        self.thresholdDict = thresholdDictionary


def _phi(listInput: np.ndarray):
    phi0 = listInput[0]
    theta0 = listInput[1]
    xi = listInput[2]
    if xi == 0:
        return xi
    else:
        return ((-2/xi)*phi0) - theta0**cf.Var.n


def _theta(listInput: np.ndarray):
    return listInput[0]


def _testLaneEmden() -> None:
    """
    Tests lane emden decoupled differential equations.
    Uses Euler for now.
    """

    # Create differential equations
    diffEq1 = DifferentialEquation(["phi", "theta", "xi"], "phi", _phi, 0, 2)
    diffEq2 = DifferentialEquation(["phi"], "theta", _theta, 1, 0)

    # Create a differential equation system.
    differentialSystem = DifferentialEquationSystem([diffEq1, diffEq2])

    # Use it to create a differential equation solver.
    differentialSolver = DifferentialSolver(differentialSystem, 0.001, 35)

    # Add a threshold for the variable theta.
    differentialSolver.addThreshold({"theta": 1e-5})

    # Iterate through all the n's and plot them
    for n in np.arange(0, 5.5, 0.5):
        cf.Var.n = n
        differentialSolver.euler()
        print(n)
        plt.plot(differentialSolver.varDict.get("xi"), differentialSolver.varDict.get("theta"))
    plt.axhline(0, color='gray', linestyle='dashed')
    plt.grid()
    plt.show()


def _RHO(inputList: np.ndarray):
    P = inputList
    return (P / cf.Var.K) ** (cf.Var.n / (cf.Var.n + 1))


def _P(inputList: np.ndarray):
    rho = inputList
    return cf.Var.K * (rho ** (1 + 1 / cf.Var.n))


def _CreateSpline():
    cubicSplineData = np.loadtxt("EoS_Spline_Data.csv", delimiter='\t').transpose()
    global _RHO
    _RHO = CubicSpline(cubicSplineData[1], cubicSplineData[0], extrapolate=True)
    global _P
    _P = CubicSpline(cubicSplineData[0], cubicSplineData[1], extrapolate=True)


def _CreateNeutronStarSpline():
    cubicSplineSlyData = np.loadtxt("wrishikData/Neutron-Star-Structure-master/SLy.txt").transpose()
    rhoData = cubicSplineSlyData[2]
    pressureData = cubicSplineSlyData[3]
    global _RHO
    _RHO = CubicSpline(pressureData, rhoData, extrapolate=True)
    global _P
    _P = CubicSpline(rhoData, pressureData, extrapolate=True)


def _M(inputList: np.ndarray):
    P = inputList[0]
    r = inputList[1]
    return 4 * np.pi * r ** 2 * _RHO(P)


def _TOV(inputList: np.ndarray):
    m = inputList[0]
    P = inputList[1]
    r = inputList[2]
    return - ((cf.G.whatUnitHuh * m * _RHO(P)) / r ** 2) \
         * (1 + P / (_RHO(P) * cf.C.cmtrPsec ** 2)) \
         * (1 + (4 * np.pi * r ** 3 * P) / (m * cf.C.cmtrPsec ** 2)) \
         * (1 - (2 * cf.G.whatUnitHuh * m) / (r * cf.C.cmtrPsec ** 2)) ** (-1)


def _HYDRO(inputList: np.ndarray):
    m = inputList[0]
    P = inputList[1]
    r = inputList[2]
    return - ((cf.G.whatUnitHuh * m * _RHO(P)) / r ** 2)


def _testTOV(rhoMin=1e6, rhoMax=1e14, rhoH=1e5, rhoNum=30, color='r', marker='v'):

    for i, rho in enumerate(np.geomspace(rhoMin, rhoMax, rhoNum)):
        diffM = DifferentialEquation(["p", "r"], "m", _M, cf.Var.m, 1)
        diffP = DifferentialEquation(["m", "p", "r"], "p", _TOV, _P(rho), 2)
        diffEqs = DifferentialEquationSystem([diffM, diffP])
        diffS = DifferentialSolver(diffEqs, rhoH, stopTime=2e11)
        rMod = diffS.varDict.get("r")
        rMod[0] = 1
        diffS.varDict.update({"r": rMod})
        diffS.addThreshold({"p": 1e-5})
        diffS.euler()
        if len(diffS.varDict.get("r")) != 0 and len(diffS.varDict.get("m")) != 0:
            r = diffS.varDict.get("r")[-1]/1e5
            m = diffS.varDict.get("m")[-1]/2e33
            print("Calculation {0} with rho = {1}\n"
                   "gave rise to r = {2} and m = {3}\n".format(i, rho, r, m))
            plt.scatter(r, m, color=color, marker=marker)
        else:
            print("Calculation {0} at rho = {1} skipped!".format(i, rho))
    # plt.show()


def _testHYDRO(rhoMin=1e6, rhoMax=1e14, rhoH=1e5, rhoNum=30, color='r', marker='v'):

    for i, rho in enumerate(np.geomspace(rhoMin, rhoMax, rhoNum)):
        diffM = DifferentialEquation(["p", "r"], "m", _M, cf.Var.m, 1)
        diffP = DifferentialEquation(["m", "p", "r"], "p", _HYDRO, _P(rho), 2)
        diffEqs = DifferentialEquationSystem([diffM, diffP])
        diffS = DifferentialSolver(diffEqs, rhoH, stopTime=2e11)
        rMod = diffS.varDict.get("r")
        rMod[0] = 1
        diffS.varDict.update({"r": rMod})
        diffS.addThreshold({"p": 1e-5})
        diffS.euler()
        if len(diffS.varDict.get("r")) != 0 and len(diffS.varDict.get("m")) != 0:
            r = diffS.varDict.get("r")[-1]/1e5
            m = diffS.varDict.get("m")[-1]/2e33
            print("Calculation {0} with rho = {1}\n"
                   "gave rise to r = {2} and m = {3}\n".format(i, rho, r, m))
            plt.scatter(r, m, color=color, marker=marker)
        else:
            print("Calculation {0} at rho = {1} skipped!".format(i, rho))
    # plt.show()


def _testN(rhoH=1e4):
    for i, n in enumerate(np.arange(0.5, 5.5, 0.5)):
        cf.Var.n = n
        diffM = DifferentialEquation(["p", "r"], "m", _M, cf.Var.m, 1)
        diffP = DifferentialEquation(["m", "p", "r"], "p", _HYDRO, _P(np.asarray([cf.Var.rho])), 2)
        diffEqs = DifferentialEquationSystem([diffM, diffP])
        diffS = DifferentialSolver(diffEqs, rhoH, stopTime=2e10)
        rMod = diffS.varDict.get("r")
        rMod[0] = 1
        diffS.varDict.update({"r": rMod})
        diffS.addThreshold({"p": 1e-5})
        diffS.euler()
        r = diffS.varDict.get("r")
        P = diffS.varDict.get("p")
        rho = _RHO(np.asarray(np.asarray(P)))
        plt.plot(r/np.max(r), rho/np.max(rho), label="n = {0}".format(cf.Var.n))
        print("Calculation {0} at n={1}".format(i, cf.Var.n))



# plt.show()

# # Polytropic
# _testTOV(rhoMin=1e7, rhoMax=1e12)
# _testHYDRO(rhoMin=1e7, rhoMax=1e12, color="b", marker="o")
# plt.grid()
# plt.show()
#
#
# # White Dwarfs
# _CreateSpline()
# _testTOV()
# _testHYDRO(color="b", marker="o")
# plt.grid()
# plt.show()
#
#
# # Neutron stars
# _CreateNeutronStarSpline()
# _testTOV(rhoMin=2e14, rhoMax=7e15, rhoH=4e3, rhoNum=30)
# _testHYDRO(rhoMin=2e14, rhoMax=1e15, rhoH=4e3, rhoNum=30, color="b", marker="o")
# plt.grid()
# plt.show()

_testN()
plt.grid()
plt.legend()
plt.show()