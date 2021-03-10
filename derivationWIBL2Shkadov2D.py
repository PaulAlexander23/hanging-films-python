#!/usr/bin/python3

from sympy import *

def main():
    wibl2TwoDimensionalFallingLiquidFilms()


def wibl2TwoDimensionalFallingLiquidFilms():
    x, y, t = symbols("x y t")
    h = Function("h")(x, t)
    q = Function("q")(x, t)
    r = Function("r")(x, t)
    s = Function("s")(x, t)
    delta, eta, zeta = symbols("delta eta zeta")
    epsilon = symbols("epsilon")

    g0 = (y/h) - (y/h)**2/Integer(2)
    g1 = (y/h) - Integer(17)/Integer(6)*(y/h)**2 + Integer(7)/Integer(3)*(y/h)**3 - Integer(7)/Integer(12)*(y/h)**4
    g2 = ((y/h) - Integer(13)/Integer(2)*(y/h)**2 + Integer(57)/Integer(4)*(y/h)**3 - Integer(111)/Integer(8)*(y/h)**4
            + Integer(99)/Integer(16)*(y/h)**5 - Integer(33)/Integer(32)*(y/h)**6)
    g = [g0, g1, g2]

    b0 = Integer(3) / h * (q - epsilon * r - epsilon * s)
    b1 = Integer(45) / h * epsilon * r
    b2 = Integer(210) / h * epsilon * s

    u = b0 * g0 + b1 * g1 + b2 * g2

    # Using the continuity equation
    vy = - diff(u, x)
    v = integrate(vy, y) # Using noslip to fix the integration constant


    # From falling liquid films
    xMomentumEquation = (
            epsilon * delta * (diff(u, t) + u * diff(u, x) + v * diff(u, y)) - Integer(2) * epsilon**2 * eta * diff(u, x, 2)
            - Integer(1) - epsilon**2 * eta * diff(diff(u, x).subs(y, h), x) + epsilon * zeta * diff(h, x) - diff(h, x, 3) 
            ).subs(diff(h, t), - diff(q, x))

    # Treat the uyy term as special
    uyy = [epsilon**2 * eta * ((4 * diff(h, x) * diff(u, x) - diff(v, x)) * gj).subs(y, h)
            + integrate(u * diff(gj, y, 2), (y, 0, h)) for gj in g]

    eqn0 = integrate(g0 * xMomentumEquation, (y, 0, h)) - uyy[0]
    eqn1 = integrate(g1 * xMomentumEquation, (y, 0, h)) - uyy[1]
    eqn2 = integrate(g2 * xMomentumEquation, (y, 0, h)) - uyy[2]

    solution = list(linsolve([eqn0, eqn1, eqn2], (diff(q, t), diff(r, t), diff(s, t))))

    #sol = [series(solution[0][0].expand(), epsilon, 0, 2), 
    #       series(solution[0][1].expand(), epsilon, 0, 1),
    #       series(solution[0][2].expand(), epsilon, 0, 1)]
    truncation = [2, 1, 1]
    sol = [series(solution[0][n].expand(), epsilon, 0, truncation[n]) for n in range(3)]

#    for n in range(3):
#        sol[n] = (sol[n].coeff(e, 0) + epsilon * sol[n].coeff(e, 1) + epsilon**2 * sol[n].coeff(e, 2)) + epsilon**3 * sol[n].coeff(e, 3))

    f = open("wibl2FallingLiquidFilms.tex","w+")
    for n in range(3):
        f.write(latex(sol[n]))
        f.write("\n")
    f.close()

    f = open("wibl2FallingLiquidFilms.txt","w+")
    for n in range(3):
        f.write(str(sol[n]))
        f.write("\n")
    f.close()

    for n in range(3):
        pprint(sol[n])
        print()

if __name__=="__main__":
    main()
