
class StringCompareError(Exception):

    def __init__(self):
        super().__init__("The to be tested substring must be smaller than the string.")


class SubStringNotInStringError(Exception):

    def __init__(self):
        super().__init__("Cannot remove substring from string, as it is not contained in string.")


def isStringPartOf(string1: str, string2: str) -> bool:
    if len(string2) > len(string1):
        raise StringCompareError

    for index, char2 in enumerate(string2):
        if string1[index] != char2:
            return False
    return True


def subtractStringBack(string1: str, string2: str) -> str:
    if isStringPartOf(string1, string2):
        return string1[len(string2) - 1: len(string1) - 1]
    else:
        raise SubStringNotInStringError
