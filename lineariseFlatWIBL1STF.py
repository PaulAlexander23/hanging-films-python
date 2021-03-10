#!/usr/bin/python3

from sympy import *
from sympy.parsing.sympy_parser import parse_expr

def main():
    x, y, z, t = symbols("x y z t")
    h = Function("h")(x, z, t)
    htilde = Function("htilde")(x, z)
    hhat = symbols("hhat")
    F1 = Function("F_1")(x, z, t)
    F1tilde = Function("F1tilde")(x, z)
    F1hat = symbols("F1hat")
    F2 = Function("F_2")(x, z, t)
    F2tilde = Function("F2tilde")(x, z)
    F2hat = symbols("F2hat")
    theta, Re, C = symbols("theta Re C")
    epsilon = symbols("epsilon")
    a, b, omega = symbols("alpha beta omega")

    ht = - diff(F1, x) - diff(F2, z)

    strings = loadStrings("wibl1-flux-eqs.txt")
    Ft = list(map(parse_expr, strings))

    eqns = [ht] + Ft
    numberOfEquations = len(eqns)

    for n in range(numberOfEquations):
        eqns[n] = eqns[n].subs(epsilon, 1)

    pert = [(h, Integer(1) + epsilon*htilde),
            (F1, Integer(2)/Integer(3) + epsilon*F1tilde),
            (F2, 0 + epsilon*F2tilde)]

    for n in range(numberOfEquations):
        eqns[n] = series(eqns[n].subs(epsilon, 1).subs(pert).doit().expand(), epsilon, 0)

    for n in range(numberOfEquations):
        eqns[n] = eqns[n].coeff(epsilon, 1)

    f = open('linearised-flat-wibl1-stf-latex.tex', 'w')
    for n in range(numberOfEquations):
        f.write(latex(eqns[n]))
        f.write("\n")
    f.close()

    phi = exp(I * (a * x + b * z))
    normalMode = [(htilde, hhat * phi),
        (F1tilde, F1hat * phi),
        (F2tilde, F2hat * phi)]

    for n in range(numberOfEquations):
        eqns[n] = powsimp((eqns[n].subs(normalMode)
            / phi).doit().expand()).expand()

    f = open('linearised-flat-wibl1-stf-exp-latex.tex', 'w')
    for n in range(numberOfEquations):
        f.write(latex(eqns[n]))
        f.write("\n")
    f.close()

    f = open('linearised-flat-wibl1-stf.txt', 'w')
    for n in range(numberOfEquations):
        f.write(str(eqns[n]))
        f.write("\n")
    f.close()

    f = open('linearised-flat-wibl1-stf-matrix.txt', 'w+')
    for n in range(numberOfEquations):
        for var in [hhat, F1hat, F2hat]:
            f.write(str(eqns[n].coeff(var, 1)))
            f.write("\n")
    f.close()

    for n in range(numberOfEquations):
        pprint(eqns[n])
        print()


def loadStrings(filename):
    f = open(filename, "r")

    stringsWithNewline = f.readlines()
    strings = list(map(str.strip, stringsWithNewline))

    f.close()

    return strings


if __name__=="__main__":
    main()
