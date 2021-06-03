#!/usr/bin/python3

from sympy import *
from sympy.parsing.sympy_parser import parse_expr
from wibl2NusseltSimplified import wibl2NusseltSimplified

def main():
    x, y, z, t = symbols("x y z t")
    h = Function("h")(x, z, t)
    hbar = Function("hbar")(z)
    htilde = Function("htilde")(x, z)
    hhat = Function("hhat")(z)
    q1 = Function("q1")(x, z, t)
    q1bar = Function("q1bar")(z)
    q1tilde = Function("q1tilde")(x, z)
    q1hat = Function("q1hat")(z)
    q2 = Function("q2")(x, z, t)
    q2bar = Function("q2bar")(z)
    q2tilde = Function("q2tilde")(x, z)
    q2hat = Function("q2hat")(z)
    theta, Re, C = symbols("theta Re C")
    epsilon = symbols("epsilon")
    delta = symbols("delta")
    alpha, beta, omega = symbols("alpha beta omega")

    ht = - diff(q1, x) - diff(q2, z)
    #Ft = wibl2NusseltSimplified(x, y, z, t, h, q1, q2, theta, Re, C)
    strings = loadStrings("wibl2STF.txt")
    Ft = list(map(parse_expr, strings))

    eqns = [ht, Ft[0], Ft[1]]

    numberOfEquations = len(eqns)

    for n in range(numberOfEquations):
        eqns[n] = eqns[n].subs([(h, hbar + delta*htilde), (q1, q1bar + delta*q1tilde), (q2, q2bar + delta*q2tilde)]).doit().expand()

    for n in range(numberOfEquations):
        eqns[n] = series(eqns[n], delta, 0)

    for n in range(numberOfEquations):
        eqns[n] = eqns[n].coeff(delta, 1)

    f = open('linearised-rivulet-wibl2-stf-latex.tex', 'w+')
    for n in range(numberOfEquations):
        f.write(latex(eqns[n]))
        f.write("\n")
    f.close

    for n in range(numberOfEquations):
        eqns[n] = powsimp((eqns[n].subs(
        [(htilde, hhat * exp(I * alpha * x + I * beta * z)),
        (q1tilde, q1hat * exp(I * alpha * x + I * beta * z)),
        (q2tilde, q2hat * exp(I * alpha * x + I * beta * z))])
            / exp(I * alpha * x + I * beta * z)).doit()).expand()

    f = open('linearised-rivulet-wibl2-stf-exp-latex.tex', 'w+')
    for n in range(numberOfEquations):
        f.write(latex(eqns[n]))
        f.write("\n")
    f.close

    f = open('linearised-rivulet-wibl2-stf.txt', 'w+')
    for n in range(numberOfEquations):
        f.write(str(eqns[n]))
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


if __name__=="__main__":
    main()
