class InvalidNumericalMethod(Exception):

    def __init__(self, wrongMethod: str):
        super().__init__("{0} is not a valid numerical method.".format(wrongMethod))


class InvalidPressureMethod(Exception):

    def __init__(self, wrongMethod):
        super().__init__("{0} is not a valid pressure correction method".format(wrongMethod))


class NoAnalyticalSolution(Exception):

    def __init__(self, wrongN):
        super().__init__("There is no valid solution for n = {o}".format(wrongN))