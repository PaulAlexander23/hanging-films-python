

from sympy import *
from saveAndLoadEquation import *

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
    delta = symbols("delta")
    alpha, beta, omega = symbols("alpha beta omega")

    strings = loadStrings("benney.txt")
    ht = list(map(parse_expr, strings))

    ht = series(ht[0].subs(h, hbar + delta*htilde).doit(), delta, 0)

    httilde = ht.coeff(delta, 1)

    f = open('linearised-rivulet-benney-latex.tex', 'w')
    f.write(latex(httilde))
    f.close

    hthat = powsimp((httilde.subs(htilde, hhat * exp(I * alpha * x + I * beta * z))
            / exp(I * alpha * x + I * beta * z)).doit().expand()).expand()

    f = open('linearised-rivulet-benney-exp-latex.tex', 'w')
    f.write(latex(hthat))
    f.close

    f = open('linearised-rivulet-benney.txt', 'w')
    f.write(str(hthat))
    f.close

    pprint(hthat)
    print()


if __name__=="__main__":
    main()
