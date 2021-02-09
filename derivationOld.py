#!/usr/bin/python3

from sympy import *

def main():
    wibl2TwoDimensionalFallingLiquidFilmsSimplified()


def wibl2TwoDimensionalFallingLiquidFilmsSimplified():
    x, y, t = symbols("x y t")
    h = Function("h")(x, t)
    q = Function("q")(x, t)
    r = Function("r")(x, t)
    s = Function("s")(x, t)
    theta, Re, C = setupParameters()
    e = symbols("epsilon")

    g0 = (y/h) - (y/h)**2/Integer(2)
    g1 = (y/h) - Integer(17)/Integer(6)*(y/h)**2 + Integer(7)/Integer(3)*(y/h)**3 - Integer(7)/Integer(12)*(y/h)**4
    g2 = ((y/h) - Integer(13)/Integer(2)*(y/h)**2 + Integer(57)/Integer(4)*(y/h)**3 - Integer(111)/Integer(8)*(y/h)**4
            + Integer(99)/Integer(16)*(y/h)**5 - Integer(33)/Integer(32)*(y/h)**6)

    b0 = Integer(3) / h * (q - e**3 * r - e**3 * s)

    u = b0 * g0

    vy = - diff(u, x)
    v = integrate(vy, y)

    # These are different
    xMomentumEquation = 3 * e * (Re * (diff(u, t) + u * diff(u, x) + v * diff(u, y)) - Integer(2) * e**2 * diff(u, x, 2) - diff(u, y, 2)
            + e * cot(theta) * diff(h, x) - C**-1 * diff(h, x, 3) - Integer(1)
            - e**2 * diff(diff(u, x).subs(y, h), x)).subs(diff(h, t), - diff(q, x))

    eqn0 = integrate(g0 * xMomentumEquation, (y, 0, h))

    solution = list(solve(eqn0, diff(q, t)))

    sol = solution[0].expand()

    sol = (sol.coeff(e, 0) + e * sol.coeff(e, 1) + e**2 * sol.coeff(e, 2))
    
    pprint(sol)


def wibl2TwoDimensionalFallingLiquidFilms():
    x, y, t = symbols("x y t")
    h = Function("h")(x, t)
    q = Function("q")(x, t)
    r = Function("r")(x, t)
    s = Function("s")(x, t)
    theta, Re, C = setupParameters()
    e = symbols("epsilon")

    g0 = (y/h) - (y/h)**2/Integer(2)
    g1 = (y/h) - Integer(17)/Integer(6)*(y/h)**2 + Integer(7)/Integer(3)*(y/h)**3 - Integer(7)/Integer(12)*(y/h)**4
    g2 = ((y/h) - Integer(13)/Integer(2)*(y/h)**2 + Integer(57)/Integer(4)*(y/h)**3 - Integer(111)/Integer(8)*(y/h)**4
            + Integer(99)/Integer(16)*(y/h)**5 - Integer(33)/Integer(32)*(y/h)**6)

    b0 = Integer(3) / h * (q - e**2 * r - e**2 * s)
    b1 = Integer(45) / h * e**2 * r
    b2 = Integer(210) / h * e**2 * s

    u = b0 * g0 + b1 * g1 + b2 * g2

    vy = - diff(u, x)
    v = integrate(vy, y)

    # These are different
    xMomentumEquation = 3 * e * (Re * (diff(u, t) + u * diff(u, x) + v * diff(u, y)) - Integer(2) * e**2 * diff(u, x, 2) - diff(u, y, 2)
            + e * cot(theta) * diff(h, x) - C**-1 * diff(h, x, 3) - Integer(1)
            - e**2 * diff(diff(u, x).subs(y, h), x)).subs(diff(h, t), - diff(q, x))

    eqn0 = integrate(g0 * xMomentumEquation, (y, 0, h))
    eqn1 = integrate(g1 * xMomentumEquation, (y, 0, h))
    eqn2 = integrate(g2 * xMomentumEquation, (y, 0, h))

    solution = list(linsolve([eqn0, eqn1, eqn2], (diff(q, t), diff(r, t), diff(s, t))))

    sol = [solution[0][n].expand() for n in range(3)]

    sol[0] = (sol[0].coeff(e, 0) + e * sol[0].coeff(e, 1) + e**2 * sol[0].coeff(e, 2))
    sol[1] = sol[1].coeff(e, 0)
    sol[2] = sol[2].coeff(e, 0)
    
    for n in range(3):
        #pprint(sol[n])
        print(latex(sol[n]))
        print()




def wibl2TwoDimensionalRuyerQuilManneville2000():
    x, y, t = symbols("x y t")
    h = Function("h")(x, t)
    q = Function("q")(x, t)
    r = Function("r")(x, t)
    s = Function("s")(x, t)
    theta, Re, C = setupParameters()
    e = symbols("epsilon")

    g0 = (y/h) - (y/h)**2/Integer(2)
    g1 = (y/h) - Integer(17)/Integer(6)*(y/h)**2 + Integer(7)/Integer(3)*(y/h)**3 - Integer(7)/Integer(12)*(y/h)**4
    g2 = ((y/h) - Integer(13)/Integer(2)*(y/h)**2 + Integer(57)/Integer(4)*(y/h)**3 - Integer(111)/Integer(8)*(y/h)**4
            + Integer(99)/Integer(16)*(y/h)**5 - Integer(33)/Integer(32)*(y/h)**6)

    b0 = Integer(3) / h * (q - e**2 * r - e**2 * s)
    b1 = Integer(45) / h * e**2 * r
    b2 = Integer(210) / h * e**2 * s

    u = b0 * g0 + b1 * g1 + b2 * g2

    vy = - diff(u, x)
    v = integrate(vy, y)

    xMomentumEquation = e * (Re * (diff(u, t) + u * diff(u, x) + v * diff(u, y)) - Integer(2) * e**2 * diff(u, x, 2) - diff(u, y, 2)
            + Integer(2) * e * cot(theta) * diff(h, x) - C**-1 * diff(h, x, 3) - Integer(2)
            - e**2 * diff(diff(u, x).subs(y, h), x)).subs(diff(h, t), - diff(q, x))

    eqn0 = integrate(g0 * xMomentumEquation, (y, 0, h))
    eqn1 = integrate(g1 * xMomentumEquation, (y, 0, h))
    eqn2 = integrate(g2 * xMomentumEquation, (y, 0, h))

    #pprint(eqn0)
    #print(degree(eqn0, gen=e))
    #eqn0 = neglectHigherOrderTerms(eqn0, e, 3)
    #pprint(eqn0)
    #print(degree(eqn0, gen=e))

    solution = list(linsolve([eqn0, eqn1, eqn2], (diff(q, t), diff(r, t), diff(s, t))))

    sol = [solution[0][n].expand() for n in range(3)]
    
    sol[0] = (sol[0].coeff(e, 0) + e * sol[0].coeff(e, 1) + e**2 * sol[0].coeff(e, 2))
    sol[1] = sol[1].coeff(e, 0)
    sol[2] = sol[2].coeff(e, 0)
    
    for n in range(3):
        #pprint(sol[n])
        print(latex(sol[n]))
        print()

    #for n in range(3):
    #    pprint((solution[0][n].expand()))
    #    #print(latex(solution[0][n].expand()))
    #    print()

    #for n in range(3):
    #    #pprint((solution[0][n].expand()).coeff(e, 0) + (solution[0][n].expand()).coeff(e, 1) + (solution[0][n].expand()).coeff(e, 2))
    #    print(solution[0][n].expand())
    #    #print(latex(solution[0][n].expand()))
    #    print()

def neglectHigherOrderTerms(expression, var, maxPower):
    return sum([var**power * expression.coeff(var, power) for power in range(maxPower)])

def wibl2ThreeDimensionalRuyerQuilManneville2006():
    x, y, z, t = symbols("x y z t")
    h = Function("h")(x, z, t)
    q1 = Function("q1")(x, z, t)
    q2 = Function("q2")(x, z, t)
    r1 = Function("r1")(x, z, t)
    r2 = Function("r2")(x, z, t)
    s1 = Function("s1")(x, z, t)
    s2 = Function("s2")(x, z, t)
    theta, Re, C = setupParameters()
    e = symbols("epsilon")

    g0 = (y/h) - (y/h)**2/Integer(2)
    g1 = (y/h) - Integer(17)/Integer(6)*(y/h)**2 + Integer(7)/Integer(3)*(y/h)**3 - Integer(7)/Integer(12)*(y/h)**4
    g2 = ((y/h) - Integer(13)/Integer(2)*(y/h)**2 + Integer(57)/Integer(4)*(y/h)**3 - Integer(111)/Integer(8)*(y/h)**4
            + Integer(99)/Integer(16)*(y/h)**5 - Integer(33)/Integer(32)*(y/h)**6)

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

    xMomentumEquation = e * (Re * (diff(u, t) + u * diff(u, x) + v * diff(u, y) + w * diff(u, z))
            - Integer(2) * e**2 * diff(u, x, 2) - diff(u, y, 2) - e**2 * diff(u, z, 2) - e**2 * diff(diff(w, x), z)
            + Integer(2) * e * cot(theta) * diff(h, x) - C**-1 * e**3 * diff(diff(h, x, 2) + diff(h, z, 2), x) - Integer(2)
            - e**2 * diff(diff(u, x).subs(y, h), x).subs(diff(h, t), - diff(q1, x) - diff(q2, x)) - e**2 * diff(diff(w, z).subs(y, h), x).subs(diff(h, t), - diff(q1, x) - diff(q2, x)))

    zMomentumEquation = e * (Re * (diff(w, t) + w * diff(w, z) + v * diff(w, y) + u * diff(w, x))
            - Integer(2) * e**2 * diff(w, z, 2) - diff(w, y, 2) - e**2 * diff(w, x, 2) - e**2 * diff(diff(u, z), x)
            + Integer(2) * e * cot(theta) * diff(h, z) - C**-1 * e**3 * diff(diff(h, z, 2) + diff(h, x, 2), z) - Integer(2)
            - e**2 * diff(diff(w, z).subs(y, h), z).subs(diff(h, t), - diff(q1, x) - diff(q2, x)) - e**2 * diff(diff(u, x).subs(y, h), z).subs(diff(h, t), - diff(q1, x) - diff(q2, x)))

    eqn0 = integrate(g0 * xMomentumEquation, (y, 0, h))
    eqn1 = integrate(g1 * xMomentumEquation, (y, 0, h))
    eqn2 = integrate(g2 * xMomentumEquation, (y, 0, h))
    eqn3 = integrate(g0 * zMomentumEquation, (y, 0, h))
    eqn4 = integrate(g1 * zMomentumEquation, (y, 0, h))
    eqn5 = integrate(g2 * zMomentumEquation, (y, 0, h))

    solution = list(linsolve([eqn0, eqn1, eqn2, eqn3, eqn4, eqn5], (diff(q1, t), diff(r1, t), diff(s1, t), diff(q2, t), diff(r2, t), diff(s2, t))))

    for n in range(6):
        print((solution[0][n].expand()))
        #print(latex(solution[0][n].expand()))
        print()


def wibl2TwoDimensional():
    x, y, t = symbols("x y t")
    h = Function("h")(x, t)
    #numberOfFunctions = 3
    #numberOfConstants = 4
    numberOfFunctions = 3
    numberOfConstants = 4
    a = ([Function("a"+str(j))(x, t) for j in range(numberOfFunctions)] 
            + [Symbol("a"+str(j+numberOfFunctions)) for j in range(numberOfConstants)])
    F1 = Function("F_1")(x, y, t)
    theta, Re, C = setupParameters()

    N = numberOfFunctions + numberOfConstants

    hxx = diff(h, x, 2)

    f = [(y/h)**(j + 1) - (Integer(j + 1))/(Integer(j + 2)) * (y/h)**(j+2) for j in range(N)]
    g = [f[0], f[1] - f[2]/Integer(3), - Integer(4) * f[2] / Integer(3) + f[3] - f[4] / Integer(5)]

    u = sum([a[j] * g[j] for j in range(N)])
    u0 = sum([a[j] * f[j] for j in [1, 3, 5]])

    vy = - diff(u0, x)
    v = integrate(vy, y)
    
    xMomentumEquation = (Re * (diff(u0, t) + u0 * diff(u0, x) + v * diff(u0, y))
            - Integer(2) * cot(theta) * diff(h, x)
            - C**-1 * diff(hxx, x)
            - diff(u, y, 2)
            - Integer(2))

    coefficientDict = Poly(xMomentumEquation, y).as_dict()

    myeqns = [coefficientDict[(deg,)] * h**deg for deg in range(N)]

    for deg in range(N):
        #print(simplify(myeqns[deg]))
        print(latex(simplify(myeqns[deg])))
        print()

    aSolution = list(linsolve(myeqns[numberOfFunctions:N], (a[numberOfFunctions:N])))
   
    for n in range(numberOfConstants):
        print(simplify(aSolution[0][n]))
        #print(latex(simplify(aSolution[0][n])))
        print()

    myeqns1 = [myeqns[n].subs([(a[numberOfFunctions+m], aSolution[0][m]) for m in range(numberOfConstants)]) * h**2 for n in range(numberOfFunctions)]

    for n in range(numberOfFunctions):
        print(simplify(myeqns1[n]))
        #print(latex(simplify(myeqns1[n])))
        print()

    atSolution = list(linsolve(myeqns1[:numberOfFunctions], (a[:numberOfFunctions])))

    for n in range(numberOfFunctions):
        print(simplify(atSolution[0][n]))
        print()


def wibl2TwoDimensionalNaive():
    x, y, t = symbols("x y t")
    h = Function("h")(x, t)
    #numberOfFunctions = 3
    #numberOfConstants = 4
    numberOfFunctions = 5
    numberOfConstants = 4
    a = ([Function("a"+str(j))(x, t) for j in range(numberOfFunctions)] 
            + [Symbol("a"+str(j+numberOfFunctions)) for j in range(numberOfConstants)])
    F1 = Function("F_1")(x, y, t)
    theta, Re, C = setupParameters()

    N = numberOfFunctions + numberOfConstants

    hxx = diff(h, x, 2)

    f = [(y/h)**(j + 1) - (Integer(j + 1))/(Integer(j + 2)) * (y/h)**(j+2) for j in range(N)]

    u = sum([a[j] * f[j] for j in range(N)])
    u0 = sum([a[j] * f[j] for j in range(numberOfFunctions)])

    vy = - diff(u0, x)
    v = integrate(vy, y)
    
    xMomentumEquation = (Re * (diff(u0, t) + u0 * diff(u0, x) + v * diff(u0, y))
            - Integer(2) * cot(theta) * diff(h, x)
            - C**-1 * diff(hxx, x)
            - diff(u, y, 2)
            - Integer(2))

    coefficientDict = Poly(xMomentumEquation, y).as_dict()

    myeqns = [coefficientDict[(deg,)] * h**deg for deg in range(N)]

    for deg in range(N):
        #print(simplify(myeqns[deg]))
        print(latex(simplify(myeqns[deg])))
        print()

    aSolution = list(linsolve(myeqns[numberOfFunctions:N], (a[numberOfFunctions:N])))
   
    for n in range(numberOfConstants):
        print(simplify(aSolution[0][n]))
        #print(latex(simplify(aSolution[0][n])))
        print()

    myeqns1 = [myeqns[n].subs([(a[numberOfFunctions+m], aSolution[0][m]) for m in range(numberOfConstants)]) * h**2 for n in range(numberOfFunctions)]

    for n in range(numberOfFunctions):
        print(simplify(myeqns1[n]))
        #print(latex(simplify(myeqns1[n])))
        print()

    atSolution = list(linsolve(myeqns1[:numberOfFunctions], (a[:numberOfFunctions])))

    for n in range(numberOfFunctions):
        print(simplify(atSolution[0][n]))
        print()


def wibl1():
    x, y, z, t = setupIndependantVariables()
    h, a, F1 = setupDependantVariables(x, y, z, t)
    theta, Re, C = setupParameters()

    hxx = diff(h, x, 2)
    hzz = diff(h, z, 2)

    f = [(y/h)**(j + 1) - (Integer(j + 1))/(Integer(j + 2)) * (y/h)**(j+2) for j in range(5)]

    u = sum([a[j] * f[j] for j in range(5)])
    u0 = a[0] * f[0]
    w1 = -(Integer(2) * cot(theta) * diff(h, z) + C**-1 * (diff(hzz, x) + diff(hzz, z))) * f[0]

    vy = - diff(u0, x) - 0 * diff(w1, z)
    v = integrate(vy, y)
    
    xMomentumEquation = (Re * (diff(u0, t) + u0 * diff(u0, x) + v * diff(u0, y) + 0 * w1 * diff(u0, z))
            - Integer(2) * cot(theta) * diff(h, x)
            - C**-1 * (diff(hxx, x) + diff(hxx, z))
            - diff(u, y, 2)
            - Integer(2))

    coefficientDict = Poly(xMomentumEquation, y).as_dict()

    myeqns = [coefficientDict[(deg,)] * h**deg for deg in range(5)]

    for deg in range(5):
        #print(simplify(myeqns[deg]))
        print(latex(simplify(myeqns[deg])))
        print()

    aSolution = list(linsolve(myeqns[1:5], (a[1:5])))
   
    myeqn1 = coefficientDict[(0,)].subs(a[1], aSolution[0][0]) * h**2

    print(latex(expand(myeqn1)))
    print()

    myeqn2 = F1 - integrate(u.subs([(a[j+1], aSolution[0][j]) for j in range(4)]),y).subs(y, h)

    print(latex(expand(myeqn2)))
    print()

    myeqn3 = myeqn1 + Integer(3) / h * myeqn2

    print(latex(expand(myeqn3)))
    print()
    
    myeqn4 = simplify(myeqn3.subs(a[0], 3 * F1 / h)).subs(diff(h, t), -diff(F1, x))

    print(latex(expand(myeqn4)))
    print()


def setupIndependantVariables():
    x, y, z, t = symbols("x y z t")
    return x, y, z, t


def setupDependantVariables(x, y, z, t):
    h = Function("h")(x, z, t)
    # Allow derivatives of 
    a = [Function("a0")(x, z, t)] + [Symbol("a"+str(j+1)) for j in range(4)]
    #a = [Function("a"+str(j))(x, z, t) for j in range(5)]
    F1 = Function("F_1")(x, y, z, t)
    return h, a, F1


def setupParameters():
    theta, Re, C = symbols("theta Re C")
    return theta, Re, C


if __name__=="__main__":
    main()
