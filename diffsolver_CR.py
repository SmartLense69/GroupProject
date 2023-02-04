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


class InputNotSortedError(Exception):
    """
    Thrown when the list of parameters are not sorted.
    """

    def __init__(self):
        super().__init__("The input variable names list is not sorted")


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
                raise InputNotSortedError()

    def call(self, inputValues: dict):
        sortedInputValues = dict(sorted(inputValues.items()))
        sizeDictionary = len(sortedInputValues)
        if self.inputVarSize != sizeDictionary:
            raise InputVariableNumberMismatchError(sizeDictionary, self.inputVarSize)
        numpyInputValues = np.zeros(sizeDictionary, dtype=np.float128)
        for _index in range(0, sizeDictionary):
            numpyInputValues[_index] = list(sortedInputValues.values())[_index]
        return self.functionCall(numpyInputValues)

    def __init__(self, inputVariables: list[str], outputVariable: str, function: callable(np.ndarray),
                 initialCondition: float):
        self.inputVarNames: list[str] = inputVariables
        self.inputVarSize: int = len(self.inputVarNames)
        self.outputVarName: str = outputVariable
        self.functionCall: callable(np.ndarray) = function
        self.initialCondit: float = initialCondition
        self.__checkInputDuplicates()
        self.__checkInputSorting()


class DifferentialEquationSolver:

    def __sortListByIO(self):
        listNewDiff = []
        for currentDiff in self.listDiffs:
            if len(listNewDiff) == 0:
                listNewDiff.append(currentDiff)
                continue
            for index, newDiff in enumerate(self.listDiffs):
                for inputVariable in newDiff.inputVarNames:
                    if currentDiff.outputVarName is inputVariable:
                        listNewDiff.insert(index, currentDiff)
        self.listDiffs = listNewDiff.copy()
        del listNewDiff

    def __init__(self, listDifferentials: list[DifferentialEquation]):
        self.listDiffs: list[DifferentialEquation] = listDifferentials
        self.__sortListByIO()


def _a(i: np.ndarray):
    return i[0]**2


def _b(i: np.ndarray):
    return i[0]


def _c(i: np.ndarray):
    return i[0] + i[1]

def _z(i: np.ndarray):
    return np.sqrt(i[0] - i[1])


diffEq1 = DifferentialEquation(["x", "y"], "z", _z, 1)
diffEq2 = DifferentialEquation(["x"], "a", _a, 2)
diffEq3 = DifferentialEquation(["a", "b"], "c", _c, 0)
diffEq4 = DifferentialEquation(["y"], "b", _b, 0)

differentialSolver = DifferentialEquationSolver([diffEq1, diffEq2, diffEq3, diffEq4])
