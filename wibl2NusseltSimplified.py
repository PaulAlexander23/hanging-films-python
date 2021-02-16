#!/usr/bin/python3

from sympy import *

def main():
    x, y, z, t = symbols("x y z t")
    h = Function("h")(x, z, t)
    F1 = Function("F1")(x, z, t)
    F2 = Function("F2")(x, z, t)
    theta, Re, C = symbols("theta Re C")

    eqns = wibl2NusseltSimplified(x, y, z, t, h, F1, F2, theta, Re, C)
    for n in range(2):
        print(latex(eqns[n]))
        print()

def wibl2NusseltSimplified(x, y, z, t, h, F1, F2, theta, Re, C):

    laph = diff(h, x, 2) + diff(h, z, 2)

    eqn1 = (-3 * Re * diff(F1, t) 
            + 3 * Re * ( Integer(9)/Integer(7) * F1**2 / h**2 * diff(h, x) 
            - Integer(17)/Integer(7) * F1 / h * diff(F1, x) ) 
            + ( Integer(5)/Integer(6) * h 
            - Integer(5)/Integer(2) * F1 / h**2 
            + 3 * Re * ( 
            - Integer(8)/Integer(7) * (F1 * diff(F2, z)) / h 
            - Integer(9)/Integer(7) * F2 * diff(F1, z) / h 
            + Integer(9)/Integer(7) * F1 * F2 * diff(h, z) / h**2) 
            + (Integer(4) * F1 * (diff(h, x))**2/h**2 
            - Integer(9)/Integer(2) *diff(F1, x) * diff(h, x)/h 
            - Integer(6) * F1 * diff(h, x, 2)/h 
            + Integer(9)/Integer(2) * diff(F1, x, 2) 
            + Integer(13)/Integer(4) *F2 * diff(h, x)* diff(h, z)/h**2 
            - diff(F1, z)*diff(h, z)/h 
            - Integer(43)/Integer(16) * diff(F2, x) * diff(h, z)/h 
            - Integer(13)/Integer(16) * diff(F2, z) * diff(h, x)/h 
            + Integer(3)/Integer(4) * F1 * (diff(h, z))**2/h**2 
            - Integer(23)/Integer(16) * F1 * diff(h, z, 2)/h 
            - Integer(73)/Integer(16) * F2 * diff(h, x, z)/h 
            + diff(F1, z, 2) 
            + Integer(7)/Integer(2) * diff(F2, x, z)) 
            - Integer(5)/Integer(6) * cot(theta) * h * diff(h, x) 
            + Integer(5)/Integer(6) * h * diff(laph, x) / C)) / 3 / Re

    eqn2 = (- 3 * Re * diff(F2, t) 
            + 3 * Re * ( Integer(9)/Integer(7) * F2**2/h**2 * diff(h, z) 
            - Integer(17)/Integer(7) * F2/h * diff(F2, z) ) 
            - Integer(5)/Integer(2) * F2/h**2 
            + 3 * Re *( 
            -Integer(8)/Integer(7) * F2 *diff(F1, x)/h 
            - Integer(9)/Integer(7) * F1 *diff(F2, x)/h 
            + Integer(9)/Integer(7) * F1 *F2 *diff(h, x)/h**2)
            + (4 * F2 *(diff(h, z))**2/h**2 
            - Integer(9)/Integer(2) *diff(F2, z) * diff(h, z)/h 
            - Integer(6) * F2 * diff(h, z, 2)/h 
            + Integer(9)/Integer(2) *diff(F2, z, 2) 
            + Integer(13)/Integer(4) * F1 * diff(h, x) * diff(h, z)/h**2 
            - diff(F2, x) * diff(h, x)/h
            - Integer(43)/Integer(16)* diff(F1, z) * diff(h, x)/h 
            - Integer(13)/Integer(16) * diff(F1, x) * diff(h, z)/h 
            + Integer(3)/Integer(4) * F2 *(diff(h, x))**2/h**2 
            - Integer(23)/Integer(16) * F2 * diff(h, x, 2)/h 
            - Integer(73)/Integer(16) * F1 * diff(h, x, z)/h 
            + diff(F2, x, 2) 
            + Integer(7)/Integer(2) * diff(F1, x, z))                        
            - Integer(5)/Integer(6) * cot(theta) * h * diff(h, z) 
            + Integer(5)/Integer(6) * h * diff(laph, z) / C) / 3 / Re

    eqns = [eqn1, eqn2]

    return eqns

if __name__=="__main__":
    main()
