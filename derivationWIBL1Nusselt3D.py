#!/usr/bin/python3

from sympy import *

def main():
    wibl1()

def wibl1():
    # Define variable
    x, y, z, t = symbols("x y z t")
    h = Function("h")(x, z, t)
    a = [Function("a0")(x, z, t)] + [Symbol("a"+str(j+1)) for j in range(4)]
    b = [Function("b0")(x, z, t)] + [Symbol("b"+str(j+1)) for j in range(4)]
    F1 = Function("F_1")(x, z, t)
    F2 = Function("F_2")(x, z, t)
    theta, Re, C = symbols("theta Re C")
    epsilon = symbols("epsilon")

    # Set up polynomials
    f = [(y/h)**(j + 1) - (Integer(j + 1))/(Integer(j + 2)) * (y/h)**(j+2) for j in range(5)]

    # Set up dependant variables
    u = sum([a[j] * f[j] for j in range(5)])
    w = sum([b[j] * f[j] for j in range(5)])
    u0 = a[0] * f[0]
    w0 = b[0] * f[0]

    # Set up equations
    v0y = - diff(u0, x) - diff(w0, z)
    v0 = integrate(v0y, y)
    
    hxx = diff(h, x, 2)
    hzz = diff(h, z, 2)

    xMomentumEquation = (epsilon * Re * (diff(u0, t) + u0 * diff(u0, x) + v0 * diff(u0, y) + w0 * diff(u0, z))
            + Integer(2) * epsilon * cot(theta) * diff(h, x)
            - C**-1 * (diff(hxx, x) + diff(hxx, z))
            - diff(u, y, 2)
            - Integer(2)
            ).subs(diff(h, t), -diff(F1, x) - diff(F2, z))

    zMomentumEquation = (epsilon * Re * (diff(w0, t) + u0 * diff(w0, x) + v0 * diff(w0, y) + w0 * diff(w0, z))
            + Integer(2) * epsilon * cot(theta) * diff(h, z)
            - C**-1 * (diff(hzz, x) + diff(hzz, z))
            - diff(w, y, 2)
            ).subs(diff(h, t), -diff(F1, x) - diff(F2, z))

    # Collect coefficients at each degree of y
    coefficientDictX = Poly(xMomentumEquation, y).as_dict()
    coefficientDictZ = Poly(zMomentumEquation, y).as_dict()

    myeqnsX = [coefficientDictX[(deg,)] * h**deg for deg in range(5)]
    myeqnsZ = [coefficientDictZ[(deg,)] * h**deg for deg in range(5)]

    myeqns = list(map(simplify, myeqnsX + myeqnsZ))

    writeLatex("wibl1-system-3d.tex", myeqns)

    # Set up linear system
    unknowns = a[1:5] + b[1:5]
    linearSystem = myeqns[1:5] + myeqns[6:10]
    aSolution = list(linsolve(linearSystem, unknowns))

    # Set up flux equations
    fluxEquation = [F1 - sum([a[j] * Integer(2) / Integer((j + 2) * (j + 3)) * h for j in range(5)]),
        F2 - sum([b[j] * Integer(2) / Integer((j + 2) * (j + 3)) * h for j in range(5)])]
    
    substitutions = list(zip(unknowns, aSolution[0]))

    for n in range(2):
        fluxEquation[n] = expand(fluxEquation[n].subs(substitutions).doit())

    printList(fluxEquation)

    writeLatex("wibl1-before-subs.tex", fluxEquation)

    writeLatex("wibl1-eqn-before-subs.tex", list(map(expand, [myeqns[0].subs(substitutions), myeqns[5].subs(substitutions)])))

    equations = [myeqns[0]*h**2 + fluxEquation[0] * Integer(3)/h, myeqns[5]*h**2 + fluxEquation[1] * Integer(3)/h]
    for n in range(2):
        equations[n] = expand(equations[n].subs(substitutions).doit())
    
    printList(equations)

    fluxSubstitution = [(a[0], Integer(3) * F1 / h), (b[0], Integer(3) * F2 / h)]
    
    for n in range(2):
        equations[n] = expand(equations[n].subs(fluxSubstitution).doit())

    for n in range(2):
        equations[n] = expand(equations[n].subs(diff(h, t), -diff(F1, x) - diff(F2, z)).doit())

    printList(equations)

    FTSolution = list(linsolve(equations, (diff(F1, t), diff(F2, t))))

    FTSolution = list(map(expand, FTSolution[0]))

    writeLatex("wibl1-flux-eqs.tex", FTSolution)

    writeStr("wibl1-flux-eqs.txt", FTSolution)

    printList(FTSolution)


def calculateResiduals(fj, xMomentumEquation, zMomentumEquation):
    linearSystemX = [integrate(fj * xMomentumEquation, (y, 0, h)) for fj in f]
    linearSystemZ = [integrate(fj * zMomentumEquation, (y, 0, h)) for fj in f]
    linearSystem = linearSystemX + linearSystemZ
    return linearSystem


def printList(varList):
    for var in varList:
        pprint(var)
        print()


def writeLatex(filename, varList):
    f = open(filename, "w+")
    for var in varList:
        f.write(latex(var))
        f.write("\n")
    print()


def writeStr(filename, varList):
    f = open(filename, "w+")
    for var in varList:
        f.write(str(var))
        f.write("\n")
    print()


if __name__=="__main__":
    main()
