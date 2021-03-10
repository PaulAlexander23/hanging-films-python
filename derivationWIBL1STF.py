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

    F1 = Function("F1")(x, z, t)
    F2 = Function("F2")(x, z, t)

    theta, Re, C = symbols("theta Re C")
    epsilon = symbols("epsilon")

    sol = list(map(parse_expr, strings))

    F1eqn = series(sol[0].subs(F2, epsilon * F2), epsilon, 0, 1).removeO()

    F2 = 

    F2sol = solve(h**2 * F2eqn, F2)

    pprint(F2sol)
    pprint(F2eqn)

    sol = F1eqn.subs(F2, F2sol[0])

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

    for n in range(len(sol)):
        pprint(sol[n])
        print()


if __name__ == "__main__":
    main()
