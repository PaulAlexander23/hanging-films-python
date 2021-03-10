#!/usr/local/bin/python3

from sympy import *
from sympy.parsing.sympy_parser import parse_expr


def main():
    strings = loadStrings("wibl2threedimensionalSimplified.txt")

    smallTransverseFlow(strings)
    

def loadStrings(filename):
    f = open(filename, "r")

    stringsWithNewline = f.readlines()
    strings = list(map(str.strip, stringsWithNewline))

    f.close()

    return strings


def smallTransverseFlow(strings):
    x, y, z, t = symbols("x y z t")
    h = Function("h")(x, z, t)

    q1 = Function("q1")(x, z, t)
    q2 = Function("q2")(x, z, t)

    theta, Re, C = symbols("theta Re C")
    epsilon = symbols("epsilon")

    sol = list(map(parse_expr, strings))

    sol = [series(sol[n].subs(q2, epsilon * q2), epsilon, 0, 2).removeO() for n in range(2)]

    

    f = open("wibl2threedimensionalSimplifiedSTF.tex","w+")
    for n in range(2):
        f.write(latex(sol[n]))
        f.write("\n")
    f.close()

    f = open("wibl2threedimensionalSimplifiedSTF.txt","w+")
    for n in range(2):
        f.write(str(sol[n]))
        f.write("\n")
    f.close()

    for n in range(len(sol)):
        pprint(sol[n])
        print()
    

if __name__ == "__main__":
    main()
