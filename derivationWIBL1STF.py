#!/usr/local/bin/python3

from sympy import *
from sympy.parsing.sympy_parser import parse_expr


def main():
    strings = loadStrings("wibl1-flux-eqs.txt")

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

    F1 = Function("F_1")(x, z, t)
    F2 = Function("F_2")(x, z, t)

    theta, Re, C = symbols("theta Re C")
    epsilon = symbols("epsilon")

    sol = list(map(parse_expr, strings))

    F1eqn = series(sol[0].subs(F2, epsilon * F2), epsilon, 0, 1).removeO()
    F2eqn = series(sol[1].subs(F2, epsilon * F2), epsilon, 0, 1).removeO()

    F2exact = solve(F2eqn, F2)

    sol = [F1eqn, F2exact[0]]

    f = open("wibl1STF.tex","w+")
    for n in range(2):
        f.write(latex(sol[n]))
        f.write("\n")
    f.close()

    f = open("wibl1STF.txt","w+")
    for n in range(2):
        f.write(str(sol[n]))
        f.write("\n")
    f.close()

    pprint(sol)
    print()


if __name__ == "__main__":
    main()
