#!/usr/bin/env python3
import numpy as np
from py3nj import wigner3j, wigner6j
from scipy.special import gamma, assoc_laguerre
from scipy import integrate
from Orbits import ElectronOrbit
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

def thj(j1, j2, j3, m1, m2, m3):
    """
    3-j symbol
    ( j1 j2 j3 )
    ( m1 m2 m3 )
    """
    return wigner3j(j1,j2,j3,m1,m2,m3)
def sjs(j1, j2, j3, j4, j5, j6):
    """
    6-j symbol
    { j1 j2 j3 }
    { j4 j5 j6 }
    """
    return wigner6j(j1,j2,j3,j4,j5,j6)

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
    print(thj(1,0,1,1,0,-1))
    print(sjs(0,0,0,0,0,0))
    print(coulomb_F_laguerre_wave_function(0,0, 0,0, 0,0, 0,0, 0,1))
