#!/usr/local/bin/python3

from sympy import *
from saveAndLoadEquation import *

def main():
    x, y, z, t = symbols("x y z t")
    h = Function("h")(x, z, t)
    htilde = Function("htilde")(x, z)
    F1 = Function("F1")(x, z, t)
    F2 = Function("F2")(x, z, t)
    P = Function("P")(x, z, t)
    theta, Re, C = symbols("theta Re C")
    epsilon = symbols("epsilon")
    delta = symbols("delta")
    a, b, omega = symbols("alpha beta omega")

    P = Integer(2) * h * cot(theta) - epsilon**2 * C**-1 * (diff(h, x, 2) + diff(h, z, 2))
    F1 = Integer(2) / Integer(3) * h**3 - epsilon * h**3 * diff(P, x) / Integer(3) + Integer(8) / Integer(15) * Re * epsilon * h**6 * diff(h, x)
    F2 = - epsilon * h**3 * diff(P, z) / Integer(3)
    ht = - diff(F1, x) - diff(F2, z)

    saveText("benney.txt",ht)

    strings = loadStrings("benney.txt")
    ht = list(map(parse_expr, strings))

    ht = series(ht[0].subs(h, 1 + delta*htilde).doit(), delta, 0)

    httilde = ht.coeff(delta, 1)

    hthat = powsimp((httilde.subs(htilde, exp(I * (a * x + b * z)))
            / exp(I * (a * x + b * z))).doit().expand()).expand()

    f = open('linearised-flat-benney-latex.tex', 'w')
    f.write(latex(hthat))
    f.close

    f = open('linearised-flat-benney.txt', 'w')
    f.write(str(hthat))
    f.close

    pprint(hthat)
    print()


def loadStrings(filename):
    f = open(filename, "r")

    stringsWithNewline = f.readlines()
    strings = list(map(str.strip, stringsWithNewline))

    f.close()

    return strings


if __name__=="__main__":
    main()
