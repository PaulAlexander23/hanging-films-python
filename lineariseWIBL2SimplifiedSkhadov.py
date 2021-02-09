#!/usr/bin/python3

from sympy import *
from wibl2SimplifiedSkhadov import wibl2SimplifiedSkhadov

def main():
    x, y, z, t = symbols("x y z t")
    h = Function("h")(x, z, t)
    hbar = Function("hbar")(z)
    htilde = Function("htilde")(x, z)
    hhat = Function("hhat")(z)
    F1 = Function("F1")(x, z, t)
    F1bar = Function("F1bar")(z)
    F1tilde = Function("F1tilde")(x, z)
    F1hat = Function("F1hat")(z)
    F2 = Function("F2")(x, z, t)
    F2bar = Function("F2bar")(z)
    F2tilde = Function("F2tilde")(x, z)
    F2hat = Function("F2hat")(z)
    delta, eta, zeta = symbols("delta eta zeta")
    e = symbols("epsilon")
    a, b, omega = symbols("alpha beta omega")

    ht = - diff(F1, x) - diff(F2, z)
    Ft = wibl2SimplifiedSkhadov(x, y, z, t, h, F1, F2, delta, eta, zeta)

    eqns = [ht, Ft[0], Ft[1]]

    numberOfEquations = len(eqns)

    for n in range(numberOfEquations):
        eqns[n] = eqns[n].subs([(h, hbar + e*htilde), (F1, F1bar + e*F1tilde), (F2, F2bar + e*F2tilde)]).doit().expand()

    for n in range(numberOfEquations):
        eqns[n] = eqns[n].coeff(e, 1)

    for n in range(numberOfEquations):
        eqns[n] = powsimp((eqns[n].subs(
        [(htilde, hhat * exp(I * (a * x - omega * t))),
        (F1tilde, F1hat * exp(I * (a * x - omega * t))),
        (F2tilde, F2hat * exp(I * (a * x - omega * t)))])
            / exp(I * (a * x - omega * t))).doit().expand()).expand()

    f = open('linearised-wibl2-latex.tex', 'w+')
    for n in range(numberOfEquations):
        f.write(latex(eqns[n]))
    f.close

    f = open('linearised-wibl2.txt', 'w+')
    for n in range(numberOfEquations):
        f.write(str(eqns[n]))
    f.close

    for n in range(numberOfEquations):
        pprint(eqns[n])


if __name__=="__main__":
    main()
