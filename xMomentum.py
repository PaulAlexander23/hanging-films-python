#!/usr/bin/python3

from sympy import *

def xMomentumEquationNusseltToShkadov():
    x, y, t = symbols("x y t")
    h = Function("h")(x, t)
    q = Function("q")(x, t)
    r = Function("r")(x, t)
    s = Function("s")(x, t)
    Ct, Re, C = symbols("Ct Re C")
    kappa, delta, eta, zeta = symbols("kappa delta eta zeta")
    e = symbols("epsilon")

    shkadovScaling = [(x, kappa*x), (t, kappa*t), (v, v/kappa),
            (Re, delta * kappa /3), (Ct, zeta * kappa)]

    xMomentumEquation = (Integer(3) * e * Re * (diff(u, t) + u * diff(u, x) + v * diff(u, y)) 
            - Integer(2) * e**2 * diff(u, x, 2) - diff(u, y, 2)
            + e * Ct * diff(h, x) - C**-1 * diff(h, x, 3) - Integer(1)
            - e**2 * diff(diff(u, x).subs(y, h), x)).subs(diff(h, t), - diff(q, x))

    pprint(xMomentumEquation.subs(shkadovScaling))


if __name__ == "__main__":
    xMomentumEquationNusseltToShkadov()
