"""
    In this module, all exceptions specific to this application are defined
"""


class InvalidNumericalMethod(Exception):
    """
        If the user defines a numerical method, that are not rk2, rk4, rkf or euler.
    """
    def __init__(self, wrongMethod: str):
        super().__init__("{0} is not a valid numerical method.".format(wrongMethod))


class InvalidPressureMethod(Exception):
    """
        If the user defines a pressure function, which is not TOV or Hydrostatic function.
    """

    def __init__(self, wrongMethod):
        super().__init__("{0} is not a valid pressure correction method".format(wrongMethod))


class NoAnalyticalSolution(Exception):
    """
        If there is no analytical solution of the Lane-Emden equation for a specified polytropic n.
    """
    def __init__(self, wrongN):
        super().__init__("There is no valid solution for n = {0}".format(wrongN))


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
                        "{0} inputs but {1} were given" \
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
