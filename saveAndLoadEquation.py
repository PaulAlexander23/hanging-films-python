#!/usr/local/bin/python3

from sympy import *

def loadStrings(filename):
    f = open(filename, "r")

    stringsWithNewline = f.readlines()
    strings = list(map(str.strip, stringsWithNewline))

    f.close()

    return strings


def saveLatex(filename, equation):
    f = open(filename, 'w')
    if isinstance(equation, list):
        for n in range(len(equation)):
            f.write(latex(equation[n]) + "\n")
    else:
        f.write(latex(equation))
    f.close()


def saveText(filename, equation):
    f = open(filename, 'w')
    if isinstance(equation, list):
        for n in range(len(equation)):
            f.write(str(equation[n]) + "\n")
    else:
        f.write(str(equation))
    f.close()
