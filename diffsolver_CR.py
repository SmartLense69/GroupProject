import numpy as np


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
        Checks if the give differential equations would be computable,
        given in the order they have been passed to the class constructor.
        :raises SortError:
            When a differential equation, which input is computed/output by another differential
            equation after the current differential equation. Example: a = func(b, c), followed by b = func(d)
            -> Throws SortError.
            Should have been: b = func(d) and then a = func(b,c)
        """
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

    def __unifyVariables(self):

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
        self.__checkSorting()
        self.__unifyVariables()


class DifferentialSolver:

    varDictInput: dict[str, int]
    varDict: dict[str, np.ndarray]

    def __resetVarDict(self):
        for variable in self.diffSystem.listAllVariables:
            self.varDict[variable] = np.zeros(self.stepsSize)
        for outputVariable, diff in zip(self.diffSystem.listOutput, self.diffSystem.listDiffs):
            self.varDict[outputVariable][0] = diff.initialCondition

    def __resetInputDict(self):
        self.varDictInput = {}
        for inputVariable in self.diffSystem.listInput:
            self.varDictInput[inputVariable] = 0

    def __init__(self, diffSystem: DifferentialEquationSystem, stepSize: float, stopTime: float):
        self.varDict = {}
        self.diffSystem = diffSystem
        self.steps = np.arange(0, stopTime+stepSize, stepSize)
        self.stepsSize = self.steps.size
        self.__resetVarDict()
        self.__resetInputDict()

    def euler(self) -> None:
        self.__resetVarDict()
        for i in range(1, self.stepsSize):
            self.__resetInputDict()
            for diff in self.diffSystem.listDiffs:
                inputList = []
                for inputVar in diff.inputVarNames:
                    inputVarCounter: int = self.varDictInput.get(inputVar)
                    if inputVarCounter == 0:
                        inputList.append(self.varDict.get(inputVar)[i - 1])
                        self.varDictInput.update({inputVar: inputVarCounter + 1})
                    else:
                        inputList.append(self.varDict.get(inputVar)[i])
                        self.varDictInput.update({inputVar: inputVarCounter + 1})
                inputList = np.asarray(inputList)
                outputArray = self.varDict.get(diff.outputVarName)
                outputArray[i] = outputArray[i - 1] + self.stepsSize * diff.functionCall(inputList)
                output = {diff.outputVarName: outputArray}
                self.varDict.update(output)


def _phi(listInput: np.ndarray):
    xi = listInput[0]
    theta0 = listInput[1]
    phi0 = listInput[2]
    if xi == 0:
        return xi
    else:
        return ((-2/xi)*phi0) - theta0


def _theta(listInput: np.ndarray):
    return listInput[0]


diffEq1 = DifferentialEquation(["phi0", "theta", "xi"], "phi", _phi, 0, 2)
diffEq2 = DifferentialEquation(["phi"], "theta", _theta, 1, 0)

differentialSystem = DifferentialEquationSystem([diffEq1, diffEq2])
differentialSolver = DifferentialSolver(differentialSystem, 1, 5)
print(differentialSolver.varDict)
differentialSolver.euler()
print(differentialSolver.varDict)
