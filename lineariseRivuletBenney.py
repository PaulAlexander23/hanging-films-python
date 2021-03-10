

from sympy import *

def main():
    x, y, z, t = symbols("x y z t")
    h = Function("h")(x, z, t)
    hbar = Function("hbar")(z)
    hhat = Function("hhat")(z)
    htilde = Function("htilde")(x, z)
    F1 = Function("F1")(x, z, t)
    F2 = Function("F2")(x, z, t)
    P = Function("P")(x, z, t)
    theta, Re, C = symbols("theta Re C")
    epsilon = symbols("epsilon")
    a, b, omega = symbols("alpha beta omega")

    P = Integer(2) * h * cot(theta) - C**-1 * (diff(h, x, 2) + diff(h, z, 2))
    F1 = Integer(2) / Integer(3) * h**3 - h**3 * diff(P, x) / Integer(3) + Integer(8) / Integer(15) * Re * h**6 * diff(h, x)
    F2 = - h**3 * diff(P, z) / Integer(3)
    ht = - diff(F1, x) - diff(F2, z)

    ht = series(ht.subs(h, hbar + epsilon*htilde).doit(), epsilon, 0)

    httilde = ht.coeff(epsilon, 1)

    hthat = powsimp((httilde.subs(htilde, hhat * exp(I * a * x))
            / exp(I * a * x)).doit().expand()).expand()

    f = open('linearised-rivulet-benney-latex.tex', 'w')
    f.write(latex(hthat))
    f.close

    f = open('linearised-rivulet-benney.txt', 'w')
    f.write(str(hthat))
    f.close

    pprint(hthat)
    print()


if __name__=="__main__":
    main()
