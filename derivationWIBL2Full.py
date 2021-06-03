#!/usr/bin/python3

from sympy import *


def main():
    wibl2()


def wibl2():
    x, y, z, t = symbols("x y z t")
    h = Function("h")(x, z, t)
    q1 = Function("q1")(x, z, t)
    q2 = Function("q2")(x, z, t)
    r1 = Function("r1")(x, z, t)
    r2 = Function("r2")(x, z, t)
    s1 = Function("s1")(x, z, t)
    s2 = Function("s2")(x, z, t)
    theta, Re, C = symbols("theta Re C")
    epsilon = symbols("epsilon")
    eta1, eta2 = symbols("eta1, eta2")

    g0 = (y/h) - (y/h)**2/Integer(2)
    g1 = (y/h) - Integer(17)/Integer(6)*(y/h)**2 + Integer(7)/Integer(3)*(y/h)**3 - Integer(7)/Integer(12)*(y/h)**4
    g2 = ((y/h) - Integer(13)/Integer(2)*(y/h)**2 + Integer(57)/Integer(4)*(y/h)**3 - Integer(111)/Integer(8)*(y/h)**4
            + Integer(99)/Integer(16)*(y/h)**5 - Integer(33)/Integer(32)*(y/h)**6)
    g = [g0, g1, g2]

    b0 = 3 / h * (q1 - r1 - s1)
    b1 = 45 / h * r1
    b2 = 210 / h * s1

    u = b0 * g0 + b1 * g1 + b2 * g2

    c0 = 3 / h * (q2 - r2 - s2)
    c1 = 45 / h * r2
    c2 = 210 / h * s2

    w = c0 * g0 + c1 * g1 + c2 * g2

    vy = - diff(u, x) - diff(w, z)
    v = integrate(vy, y)


    # Treat the uyy and wyy term as special
    uyy = [epsilon**2 * ((Integer(2) * diff(h, x) * (Integer(2) * diff(u, x) + diff(w, z)) + diff(h, z) * (diff(u, z) + diff(w, x)) - diff(v, x)) * gj).subs(y, h)
            + integrate(u * diff(gj, y, 2), (y, 0, h)) for gj in g]

    wyy = [epsilon**2 * ((Integer(2) * diff(h, z) * (diff(u, x) + Integer(2) * diff(w, z)) + diff(h, x) * (diff(u, z) + diff(w, x)) - diff(v, z)) * gj).subs(y, h)
            + integrate(u * diff(gj, y, 2), (y, 0, h)) for gj in g]

    xMomentumEquation = epsilon * (Re * (diff(u, t) + u * diff(u, x) + v * diff(u, y) + w * diff(u, z))
            - Integer(2) * epsilon**2 * diff(u, x, 2) - epsilon**2 * diff(u, z, 2) - epsilon**2 * diff(diff(w, x), z)
            + Integer(2) * epsilon * cot(theta) * diff(h, x) - C**-1 * epsilon**3 * diff(diff(h, x, 2) + diff(h, z, 2), x) - Integer(2) * eta1
            - epsilon**2 * diff(diff(u, x).subs(y, h), x).subs(diff(h, t), - diff(q1, x) - diff(q2, x)) - epsilon**2 * diff(diff(w, z).subs(y, h), x).subs(diff(h, t), - diff(q1, x) - diff(q2, z)))

    zMomentumEquation = epsilon * (Re * (diff(w, t) + w * diff(w, z) + v * diff(w, y) + u * diff(w, x))
            - Integer(2) * epsilon**2 * diff(w, z, 2) - epsilon**2 * diff(w, x, 2) - epsilon**2 * diff(diff(u, z), x)
            + Integer(2) * epsilon * cot(theta) * diff(h, z) - C**-1 * epsilon**3 * diff(diff(h, z, 2) + diff(h, x, 2), z) - Integer(2) * eta2
            - epsilon**2 * diff(diff(w, z).subs(y, h), z).subs(diff(h, t), - diff(q1, x) - diff(q2, x)) - epsilon**2 * diff(diff(u, x).subs(y, h), z).subs(diff(h, t), - diff(q1, x) - diff(q2, z)))

    eqn0 = integrate(g0 * xMomentumEquation, (y, 0, h)) - uyy[0]
    eqn1 = integrate(g1 * xMomentumEquation, (y, 0, h)) - uyy[1]
    eqn2 = integrate(g2 * xMomentumEquation, (y, 0, h)) - uyy[2]
    eqn3 = integrate(g0 * zMomentumEquation, (y, 0, h)) - wyy[0]
    eqn4 = integrate(g1 * zMomentumEquation, (y, 0, h)) - wyy[1]
    eqn5 = integrate(g2 * zMomentumEquation, (y, 0, h)) - wyy[2]

    solution = list(linsolve([eqn0, eqn1, eqn2, eqn3, eqn4, eqn5], (diff(q1, t), diff(r1, t), diff(s1, t), diff(q2, t), diff(r2, t), diff(s2, t))))

    truncation = [2, 1, 1, 2, 1, 1]
    sol = [series(solution[0][n].expand(), epsilon, 0, truncation[n]) for n in range(6)]

    sol = [expand(sol[n].subs(diff(h,t),-diff(q1,x)-diff(q2,z))) for n in range(6)]

    f = open("wibl2threedimensional.tex","w+")
    for n in range(6):
        f.write(latex(sol[n]))
        f.write("\n")
    f.close()

    f = open("wibl2threedimensional.txt","w+")
    for n in range(6):
        f.write(str(sol[n]))
        f.write("\n")
    f.close()

    for n in range(6):
        pprint(sol[n])
        print()


if __name__ == "__main__":
    main()
