
#!/usr/bin/python3

from sympy import *

def main():
    wibl1()

def wibl1():
    x, y, z, t = symbols("x y z t")
    h, a, F1 = setupDependantVariables(x, y, z, t)
    theta, Re, C = symbols("theta Re C")
    epsilon = symbols("epsilon")

    hxx = diff(h, x, 2)
    hzz = diff(h, z, 2)

    f = [(y/h)**(j + 1) - (Integer(j + 1))/(Integer(j + 2)) * (y/h)**(j+2) for j in range(5)]

    u = sum([a[j] * f[j] for j in range(5)])
    u0 = a[0] * f[0]
    w1 = -(Integer(2) * cot(theta) * diff(h, z) + C**-1 * (diff(hzz, x) + diff(hzz, z))) * f[0]

    vy = - diff(u0, x) - 0 * diff(w1, z)
    v = integrate(vy, y)
    
    xMomentumEquation = (epsilon * Re * (diff(u0, t) + u0 * diff(u0, x) + v * diff(u0, y) + 0 * w1 * diff(u0, z))
            - Integer(2) * epsilon * cot(theta) * diff(h, x)
            - C**-1 * (diff(hxx, x) + diff(hxx, z))
            - diff(u, y, 2)
            - Integer(2))

    coefficientDict = Poly(xMomentumEquation, y).as_dict()

    myeqns = [coefficientDict[(deg,)] * h**deg for deg in range(5)]

    aSolution = list(linsolve(myeqns[1:5], (a[1:5])))
   
    #myeqn1 = coefficientDict[(0,)].subs(a[1], aSolution[0][0]) * h**2
    #pprint(expand(myeqn1))
    #print()
    #myeqn2 = F1 - integrate(u.subs([(a[j+1], aSolution[0][j]) for j in range(4)]),y).subs(y, h)
    #pprint(expand(myeqn2))
    #print()
    #myeqn3 = myeqn1 + Integer(3) / h * myeqn2
    #pprint(expand(myeqn3))
    #print()
    #myeqn4 = expand(simplify(myeqn3.subs(a[0], 3 * F1 / h)).subs(diff(h, t), -diff(F1, x)))

    myeqn1 = coefficientDict[(0,)].subs(a[1], aSolution[0][0]) * h**2
    pprint(expand(myeqn1))
    print()
    myeqn2 = Integer(3) / h * (F1 - integrate(u.subs([(a[j+1], aSolution[0][j]) for j in range(4)]),y).subs(y, h))
    pprint(expand(myeqn2))
    print()
    myeqn3 = myeqn1 + myeqn2
    pprint(expand(myeqn3))
    print()
    myeqn4 = expand(simplify(myeqn3.subs(a[0], 3 * F1 / h)).subs(diff(h, t), -diff(F1, x)))

    f = open("wibl1.txt", "w+")
    f.write(str(myeqn4))
    print()

    f = open("wibl1.tex", "w+")
    f.write(latex(myeqn4))
    print()

    pprint(expand(myeqn4))
    print()


def setupDependantVariables(x, y, z, t):
    h = Function("h")(x, z, t)
    # Allow derivatives of 
    a = [Function("a0")(x, z, t)] + [Symbol("a"+str(j+1)) for j in range(4)]
    #a = [Function("a"+str(j))(x, z, t) for j in range(5)]
    F1 = Function("F_1")(x, y, z, t)
    return h, a, F1


if __name__=="__main__":
    main()
