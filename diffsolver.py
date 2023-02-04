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
        self._message = "The differential equation requires" \
                        "{0} inputs but {1} were given"\
            .format(self._lenInputVariables, self._lenInputDictionary)
        super().__init__(self._message)


class InputNotSortedError(Exception):
    """
    Thrown when the list of parameters are not sorted.
    """

    def __init__(self):
        super().__init__("The input variable names list is not sorted")


class DiffEquation:

    def _checkInputDuplicates(self):
        for _index in range(0, self._inputVarSize - 1):
            _currentInputVarName: str = self._inputVarNames[_index]
            for _compareIndex in range(_index+1, self._inputVarSize):
                if _currentInputVarName == self._inputVarNames[_compareIndex]:
                    raise InputVariableDuplicateError(_currentInputVarName)

    def _checkInputSorting(self):
        _sortedInput = self._inputVarNames.copy()
        _sortedInput.sort()
        for _sortedInputItem, _inputVarNamesItem in zip(_sortedInput, self._inputVarNames):
            if _sortedInputItem is not _inputVarNamesItem:
                raise InputNotSortedError()

    def call(self, inputValues: dict):
        _sortedInputValues = dict(sorted(inputValues.items()))
        _sizeDictionary = len(_sortedInputValues)
        if self._inputVarSize != _sizeDictionary:
            raise InputVariableNumberMismatchError(_sizeDictionary, self._inputVarSize)
        _numpyInputValues = np.zeros(_sizeDictionary, dtype=np.float128)
        for _index in range(0, _sizeDictionary):
            _numpyInputValues[_index] = list(_sortedInputValues.values())[_index]
        return self._functionCall(_numpyInputValues)

    def __init__(self, inputVariables: list[str], outputVariable: str, function: callable(np.ndarray),
                 initialCondition: float):
        self._inputVarNames: list[str] = inputVariables
        self._inputVarSize: int = len(self._inputVarNames)
        self._outputVarName: str = outputVariable
        self._functionCall: callable(np.ndarray) = function
        self._initialCondit: float = initialCondition
        self._checkInputDuplicates()
        self._checkInputSorting()


def theta(var: np.ndarray):
    x = var[0]
    y = var[1]
    z = var[2]
    return x + y + z


diffEquation = DiffEquation(["x", "y", "z"], "theta", theta, 1)
result = diffEquation.call({"x": 2, "y": 1, "z": 5})
print(result)