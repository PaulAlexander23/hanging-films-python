#!/usr/bin/python3

from sympy import *
from sympy.parsing.sympy_parser import parse_expr
from wibl2NusseltSimplified import wibl2NusseltSimplified

def main():
    x, y, z, t = symbols("x y z t")
    h = Function("h")(x, z, t)
    htilde = Function("htilde")(x, z)
    F1 = Function("q1")(x, z, t)
    F1tilde = Function("F1tilde")(x, z)
    F2 = Function("q2")(x, z, t)
    F2tilde = Function("F2tilde")(x, z)
    hhat, F1hat, F2hat = symbols("hhat F1hat F2hat")
    theta, Re, C = symbols("theta Re C")
    e = symbols("epsilon")
    a, b, omega = symbols("alpha beta omega")

    #ht = - diff(F1, x) - diff(F2, z)
    #Ft = wibl2NusseltSimplified(x, y, z, t, h, F1, F2, theta, Re, C)

    ht = - diff(F1, x) - diff(F2, z)
    strings = loadStrings("wibl2threedimensionalSimplified.txt")
    Ft = parseStrings(strings)

    eqns = [ht, Ft[0].subs(e, 1), Ft[1].subs(e, 1)]

    numberOfEquations = len(eqns)

    for n in range(numberOfEquations):
        eqns[n] = eqns[n].subs([(h, Integer(1) + e*htilde), (F1, Integer(2)/Integer(3) + e*F1tilde), (F2, e**2*F2tilde)]).doit().expand()

    for n in range(numberOfEquations):
        eqns[n] = series(eqns[n], e)

    for n in range(numberOfEquations):
        eqns[n] = eqns[n].coeff(e, 1)

    #for n in range(numberOfEquations):
    #    pprint(eqns[n])

    f = open('linearised-wibl2-flat-stf-latex.tex', 'w+')
    for n in range(numberOfEquations):
        f.write(latex(eqns[n]))
        f.write("\n")

    for n in range(numberOfEquations):
        eqns[n] = powsimp((eqns[n].subs(
        [(htilde, hhat * exp(I * (a * x + b * z - omega * t))),
        (F1tilde, F1hat * exp(I * (a * x + b * z - omega * t))),
        (F2tilde, F2hat * exp(I * (a * x + b * z - omega * t)))])
            / exp(I * (a * x + b * z - omega * t))).doit()).expand()

    for n in range(numberOfEquations):
        f.write(latex(eqns[n]))
        f.write("\n")
    f.close

    f = open('linearised-wibl2-flat-stf.txt', 'w+')
    for n in range(numberOfEquations):
        f.write(str(eqns[n]))
        f.write("\n")
    f.close

    f = open('linearised-wibl2-flat-stf-matrix.txt', 'w+')
    for n in range(numberOfEquations):
        for var in [hhat, F1hat, F2hat]:
            f.write(str(eqns[n].coeff(var, 1)))
            f.write("\n")
    f.close

    for n in range(numberOfEquations):
        pprint(eqns[n])


def loadStrings(filename):
    f = open(filename, "r")

    stringsWithNewline = f.readlines()
    strings = list(map(str.strip, stringsWithNewline))

    f.close()

    return strings


def parseStrings(strings):
    sol = list(map(parse_expr, strings))
    return sol


if __name__=="__main__":
    main()
