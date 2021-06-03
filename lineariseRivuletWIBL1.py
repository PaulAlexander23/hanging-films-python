#!/usr/bin/python3

from sympy import *
from sympy.parsing.sympy_parser import parse_expr

def main():
    x, y, z, t = symbols("x y z t")
    h = Function("h")(x, z, t)
    hbar = Function("hbar")(z)
    htilde = Function("htilde")(x, z)
    hhat = Function("hhat")(z)
    F1 = Function("F_1")(x, z, t)
    F1bar = Function("F1bar")(z)
    F1tilde = Function("F1tilde")(x, z)
    F1hat = Function("F1hat")(z)
    F2 = Function("F_2")(x, z, t)
    F2bar = Function("F2bar")(z)
    F2tilde = Function("F2tilde")(x, z)
    F2hat = Function("F2hat")(z)
    theta, Re, C = symbols("theta Re C")
    epsilon = symbols("epsilon")
    delta = symbols("delta")
    alpha, beta, omega = symbols("alpha beta omega")

    ht = - diff(F1, x) - diff(F2, z)

    strings = loadStrings("wibl1-flux-eqs.txt")
    Ft = list(map(parse_expr, strings))

    eqns = [ht] + Ft
    numberOfEquations = len(eqns)

    pert = [(h, hbar + delta*htilde), (F1, F1bar + delta*F1tilde), (F2, F2bar + delta*F2tilde)]

    for n in range(numberOfEquations):
        eqns[n] = series(eqns[n].subs(pert).doit().expand(), delta, 0)

    for n in range(numberOfEquations):
        eqns[n] = eqns[n].coeff(delta, 1)

    f = open('linearised-rivulet-wibl1-tilde-latex.tex', 'w')
    for n in range(numberOfEquations):
        f.write(latex(eqns[n]))
        f.write("\n")
    f.close()

    for n in range(numberOfEquations):
        pprint(eqns[n])
        print()
    print()

    normalMode = [(htilde, hhat * exp(I * alpha * x + I * beta * z)),
        (F1tilde, F1hat * exp(I * alpha * x + I * beta * z)),
        (F2tilde, F2hat * exp(I * alpha * x + I * beta * z))]

    for n in range(numberOfEquations):
        eqns[n] = powsimp((eqns[n].subs(normalMode)
            / exp(I * alpha * x + I * beta * z)).doit().expand()).expand()

    f = open('linearised-rivulet-wibl1-hat-latex.tex', 'w')
    for n in range(numberOfEquations):
        f.write(latex(eqns[n]))
        f.write("\n")
    f.close()

    f = open('linearised-rivulet-wibl1.txt', 'w')
    for n in range(numberOfEquations):
        f.write(str(eqns[n]))
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
