#!/usr/bin/python3

from sympy import *

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
    F2 = Function("F2")(x, z, t)
    P = Function("P")(x, z, t)
    theta, Re, C = symbols("theta Re C")
    epsilon = symbols("epsilon")
    delta = symbols("delta")
    alpha, beta, omega = symbols("alpha beta omega")

    #P = Integer(2) * h * cot(theta) - C**-1 * (diff(h, x, 2) + diff(h, z, 2))
    #F2 = - h**3 * diff(P, z) / Integer(3)
    #ht = - diff(F1, x) - diff(F2, z)
    #F1t = Integer(5) / Integer(3) / Re * h - Integer(5) / Integer(6) / Re * h * diff(P, x) - Integer(5) /Integer(2) / Re * F1 / h**2 + Integer(9) / Integer(7) * diff(h, x) * F1**2 / h**2 - Integer(17) / Integer(7) * F1 * diff(F1, x) / h

    #eqns = [ht, F1t]
    #numberOfEquations = len(eqns)

    ht = - diff(F1, x) - diff(F2, z)
    strings = loadStrings("wibl1STF.txt")
    Ft = list(map(parse_expr, strings))

    eqns = [ht.subs(F2, Ft[1]), Ft[0]]
    numberOfEquations = len(eqns)

    for n in range(numberOfEquations):
        eqns[n] = series(eqns[n].subs([(h, hbar + delta*htilde), (F1, F1bar + delta*F1tilde)]).doit().expand(), delta, 0)

    for n in range(numberOfEquations):
        eqns[n] = eqns[n].coeff(delta, 1)

    for n in range(numberOfEquations):
        pprint(eqns[n])
        print()

    f = open('linearised-rivulet-wibl1-stf-latex.tex', 'w')
    for n in range(numberOfEquations):
        f.write(latex(eqns[n]))
        f.write("\n")
    f.close()

    for n in range(numberOfEquations):
        eqns[n] = powsimp((eqns[n].subs(
        [(htilde, hhat * exp(I * alpha * x + I * beta * z)),
        (F1tilde, F1hat * exp(I * alpha * x + I * beta * z))])
            / exp(I * alpha * x + I * beta * z)).doit().expand()).expand()

    f = open('linearised-rivulet-wibl1-stf-exp-latex.tex', 'w')
    for n in range(numberOfEquations):
        f.write(latex(eqns[n]))
        f.write("\n")
    f.close()

    f = open('linearised-rivulet-wibl1-stf.txt', 'w')
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
