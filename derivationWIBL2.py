#!/usr/local/bin/python3

from sympy import *


def main():
    wibl2Simplified()
    

def wibl2Simplified():
    x, y, z, t = symbols("x y z t")
    h = Function("h")(x, z, t)

    q1 = Function("q1")(x, z, t)
    q2 = Function("q2")(x, z, t)

    theta, Re, C = symbols("theta Re C")
    epsilon = symbols("epsilon")

    g0 = (y/h) - (y/h)**2/Integer(2)
    g1 = ((y/h) - Integer(17)/Integer(6)*(y/h)**2 
            + Integer(7)/Integer(3)*(y/h)**3 - Integer(7)/Integer(12)*(y/h)**4)
    g2 = ((y/h) - Integer(13)/Integer(2)*(y/h)**2 
            + Integer(57)/Integer(4)*(y/h)**3 - Integer(111)/Integer(8)*(y/h)**4
            + Integer(99)/Integer(16)*(y/h)**5 - Integer(33)/Integer(32)*(y/h)**6)

    g = [g0, g1, g2]

    b0 = Integer(3) / h * (q1)
    c0 = Integer(3) / h * (q2)

    u = b0 * g0
    w = c0 * g0

    vy = - diff(u, x) - diff(w, z)
    v = integrate(vy, y)


    # Treat the uyy and wyy term as special
    uyy = epsilon**2 * ((Integer(2) * diff(h, x) * (Integer(2) * diff(u, x) + diff(w, z)) + diff(h, z) * (diff(u, z) + diff(w, x)) - diff(v, x)) * g0).subs(y, h) - q1 / h**2
    wyy = epsilon**2 * ((Integer(2) * diff(h, z) * (diff(u, x) + Integer(2) * diff(w, z)) + diff(h, x) * (diff(u, z) + diff(w, x)) - diff(v, z)) * g0).subs(y, h) - q2 / h**2

    xMomentumEquation = (
            epsilon * Re * (diff(u, t) + u * diff(u, x) + v * diff(u, y) + w * diff(u, z))
            - Integer(2) * epsilon**2 * diff(u, x, 2) - epsilon**2 * diff(u, z, 2) - epsilon**2 * diff(diff(w, x), z)
            + Integer(2) * epsilon * cot(theta) * diff(h, x) - C**-1 * epsilon**3 * diff(diff(h, x, 2) + diff(h, z, 2), x) - Integer(2)
            - epsilon**2 * diff(diff(u, x).subs(y, h), x) - epsilon**2 * diff(diff(w, z).subs(y, h), x)
            ).subs(diff(h, t), - diff(q1, x) - diff(q2, z))

    zMomentumEquation = (
            epsilon * Re * (diff(w, t) + w * diff(w, z) + v * diff(w, y) + u * diff(w, x))
            - Integer(2) * epsilon**2 * diff(w, z, 2) - epsilon**2 * diff(w, x, 2) - epsilon**2 * diff(diff(u, z), x)
            + Integer(2) * epsilon * cot(theta) * diff(h, z) - C**-1 * epsilon**3 * diff(diff(h, z, 2) + diff(h, x, 2), z)
            - epsilon**2 * diff(diff(w, z).subs(y, h), z) - epsilon**2 * diff(diff(u, x).subs(y, h), z)
            ).subs(diff(h, t), - diff(q1, x) - diff(q2, z))

    eqn0 = integrate(g0 * xMomentumEquation, (y, 0, h)) - uyy
    eqn3 = integrate(g0 * zMomentumEquation, (y, 0, h)) - wyy

    solution = list(linsolve([eqn0, eqn3], (diff(q1, t), diff(q2, t))))

    sol = [solution[0][n].subs(C, epsilon**3 * C).expand() for n in range(2)]

    f = open("wibl2threedimensionalSimplified.tex","w+")
    for n in range(2):
        f.write(latex(sol[n]))
        f.write("\n")
    f.close()

    f = open("wibl2threedimensionalSimplified.txt","w+")
    for n in range(2):
        f.write(str(sol[n]))
        f.write("\n")
    f.close()

    for n in range(2):
        pprint(sol[n])
        print()


if __name__ == "__main__":
    main()
