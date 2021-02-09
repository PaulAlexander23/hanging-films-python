#!/usr/bin/python3

from sympy import *

def main():
    x, y, z, t = symbols("x y z t")
    h = Function("h")(x, z, t)
    hbar = Function("hbar")(z)
    eta = Function("eta")(z)
    htilde = Function("htilde")(x, z)
    F1 = Function("F1")(x, z, t)
    F2 = Function("F2")(x, z, t)
    P = Function("P")(x, z, t)
    theta, Re, C = symbols("theta Re C")
    e = symbols("epsilon")
    a, b, omega = symbols("alpha beta omega")

    P = Integer(2) * h * cot(theta) - C**-1 * (diff(h, x, 2) + diff(h, z, 2))
    F1 = Integer(2) / Integer(3) * h**3 - h**3 * diff(P, x) / Integer(3) + Integer(8) / Integer(15) * h**6 * diff(h, x)
    F2 = - h**3 * diff(P, z) / Integer(3)
    ht = diff(F1, x) + diff(F2, z)

    ht = ht.subs(h, hbar + e*htilde).doit().expand()

    httilde = ht.coeff(e, 1)

    hthat = powsimp((httilde.subs(htilde, eta * exp(I * (a * x - omega * t))) 
            / exp(I * (a * x - omega * t))).doit().expand()).expand()

    f = open('linearised-benney-latex.tex', 'w')
    f.write(latex(hthat))
    f.close

    f = open('linearised-benney.txt', 'w')
    f.write(str(hthat))
    f.close


if __name__=="__main__":
    main()
