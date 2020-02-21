#!/usr/bin/env python3
import numpy as np
#from py3nj import wigner3j, wigner6j
from sympy.physics.wigner import wigner_3j, wigner_6j
from sympy import N
from scipy.special import gamma, assoc_laguerre, eval_gegenbauer
from scipy import integrate
from Orbits import ElectronOrbit

c = 137.035999084
g = 2.00231930436256
def triangle(J1,J2,J3):
    """
    triangular condition
    """
    b = True
    if(abs(J1-J2) <= J3 <= J1+J2): b = False
    return b

def laguerre_wave_function(x, zeta, n, l):
    """
    Laguerre function, see [A. E. McCoy and M. A. Caprio, J. Math. Phys. 57, (2016).] for details
    """
    eta = 2.0 * x / zeta
    return np.sqrt(2.0 * gamma(n+1) / (zeta * gamma(n+2*l+3)) ) * 2.0 * eta**l * np.exp(-0.5*eta) * assoc_laguerre(eta, n, 2*l+2) / zeta

def laguerre_wave_function_mom(p, zeta, n, l):
    """
    Brute force fourier transformation of Laguerre function
    """
    xi = zeta * p * 0.5
    x = 0.5 / np.sqrt( 0.25 + xi*xi )
    r = 0.0
    prefact = np.sqrt( zeta * gamma(n+1) * gamma(n+2*l+3) /np.pi ) * (0.5*zeta) * (2*x)**(2*l+3) * (2*xi)**l * gamma(l+1)
    for k in range(n+1):
        r += (-1)**k * (k+1) * (2*x)**k * eval_gegenbauer(k+1, l+1, x) / (gamma(n+1-k) * gamma(2*l+3+k) )
    return r*prefact

def thj(j1, j2, j3, m1, m2, m3):
    """
    3-j symbol
    ( j1 j2 j3 )
    ( m1 m2 m3 )
    """
    #return wigner3j(j1,j2,j3,m1,m2,m3)
    return N(wigner_3j(0.5*j1,0.5*j2,0.5*j3,0.5*m1,0.5*m2,0.5*m3))
def sjs(j1, j2, j3, j4, j5, j6):
    """
    6-j symbol
    { j1 j2 j3 }
    { j4 j5 j6 }
    """
    #return wigner6j(j1,j2,j3,j4,j5,j6)
    return N(wigner_6j(0.5*j1,0.5*j2,0.5*j3,0.5*j4,0.5*j5,0.5*j6))

def norm_laguerre_wave_function_mom(n1, n2, l):
    """
    """
    return integrate.quad(lambda x: laguerre_wave_function_mom(x, 1.0, n1, l) * laguerre_wave_function_mom(x, 1.0, n2, l) * x*x, 0, np.inf)

def rel_laguerre_wave_function_int(n1, n2, l):
    """
    """
    return integrate.quad(lambda x: laguerre_wave_function_mom(x, 1.0, n1, l) * laguerre_wave_function_mom(x, 1.0, n2, l) * x*x*x*x*x*x*0.125/(c*c), 0, np.inf)[0]


def T_laguerre_wave_function(n1, n2, l):
    """
    Kinetic term with Laguerre function basis, see [A. E. McCoy and M. A. Caprio, J. Math. Phys. 57, (2016).] for details
    """
    if(n1>n2):
        return np.sqrt( gamma(n1+1) * gamma(n2+2*l+3) / ( gamma(n2+1) * gamma(n1+2*l+3)) ) * (4*n2+4*l+6) / (2*l+3) * 0.5
    if(n1==n2):
        return (4*n1+2*l+3) / (2*l+3) * 0.5
    if(n1<n2):
        return np.sqrt( gamma(n2+1) * gamma(n1+2*l+3) / ( gamma(n1+1) * gamma(n2+2*l+3)) ) * (4*n1+4*l+6) / (2*l+3) * 0.5

def Coul_laguerre_wave_function1(n1, n2, l, zeta):
    """
    Nucleus-electron potential
    """
    return integrate.quad(lambda x: laguerre_wave_function(x, zeta, n1, l) * laguerre_wave_function(x, zeta, n2, l) * x, 0, np.inf)

def darwin_term(n1, n2, zeta):
    return laguerre_wave_function(0, zeta, n1, 0) * laguerre_wave_function(0, zeta, n2, 0)/(8*c*c)

def spin_orb(n1, n2, l, j, zeta):
    integral = integrate.quad(lambda x: laguerre_wave_function(x, zeta, n1, l) * laguerre_wave_function(x, zeta, n2, l)/(x), 1e-10, np.inf)[0]
    return g*(0.5*j*(0.5*j + 1) - l*(l+1) - 0.75)*integral/(4*c*c)

def coulomb_F_laguerre_wave_function(n1,l1,n2,l2,n3,l3,n4,l4,L, zeta):
    """
    electron-electron integral
    """
    return integrate.dblquad(lambda r1, r2: \
            laguerre_wave_function(r1,zeta,n1,l1) * laguerre_wave_function(r2,zeta,n2,l2) * \
            laguerre_wave_function(r1,zeta,n3,l3) * laguerre_wave_function(r2,zeta,n4,l4) * \
            r1*r1*r2*r2 * min(r1,r2)**L / max(r1,r2)**(L+1), 0, np.inf, lambda r1: 0, lambda r1: np.inf)

def ee_laguerre_wave_function(o1,o2,o3,o4,J,zeta):
    """
    electron-electron interaction
    """
    r = 0.0
    lmin = max( abs(o1.j-o3.j), abs(o2.j-o4.j) )//2
    lmax = min(     o1.j+o3.j ,     o2.j+o4.j  )//2
    integral = 0.0
    for l in range(lmin,lmax+1):
        if( (o1.l+o3.l+l)%2 == 1): continue
        if( (o2.l+o4.l+l)%2 == 1): continue
        integral = coulomb_F_laguerre_wave_function(o1.n, o1.l, o2.n, o2.l, o3.n, o3.l, o4.n, o4.l, l, zeta)[0]
        r += integral * sjs(o1.j, o2.j, 2*J, o4.j, o3.j, 2*l) * thj(o1.j, 2*l, o3.j, -1, 0, 1) * thj(o2.j, 2*l, o4.j, -1, 0, 1)
    return r * np.sqrt( (o1.j+1) * (o2.j+1) * (o3.j+1) * (o4.j+1) ) * (1-2*( ( (o1.j+o3.j)/2+J)%2 ) )
if(__name__=="__main__"):
    #for x in np.arange(0.0, 10.0, 1):
    #    print(x, laguerre_wave_function(x, 1.0, 0, 0))
    l = 0
    #for n1 in range(5):
    #    t2 = norm_laguerre_wave_function_mom(n1,n1,l)[0]
    #    print(n1, n1, t2)
    for n1 in range(5):
        for n2 in range(5):
            t = T_laguerre_wave_function(n1,n2,l)
            t2 = T_laguerre_wave_function_int(n1,n2,l)[0]
            print(n1, n2, t-t2)
    #for l in range(5):
    #    print( laguerre_wave_function(0.0, 1.0, 3, l) )

    #print(thj(1,0,1,1,0,-1))
    #print(sjs(0,0,0,0,0,0))
    #print(coulomb_F_laguerre_wave_function(0,0, 0,0, 0,0, 0,0, 0,1))
