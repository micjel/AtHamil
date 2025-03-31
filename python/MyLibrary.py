#!/usr/bin/env python3
import numpy as np
import math
import time
#from py3nj import wigner3j, wigner6j
from sympy.physics.wigner import wigner_3j, wigner_6j, wigner_9j
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
    return math.sqrt(2.0 * gamma(n+1) / (zeta * gamma(n+2*l+3)) ) * 2.0 * eta**l * np.exp(-0.5*eta) * assoc_laguerre(eta, n, 2*l+2) / zeta

def dlaguerre_wave_function(x, zeta, n, l):
    """
    r-derivative of Laguerre radial function
    """
    eta = 2.0 * x / zeta
    prefact = 2.0 / zeta * math.sqrt(2.0 * gamma(n+1) / (zeta * gamma(n+2*l+3)) )
    
    if l == 0:
        return prefact*(eta**l * (-1.0) * np.exp(-0.5*eta) * assoc_laguerre(eta, n, 2*l+2) + eta**l * np.exp(-0.5*eta) * assoc_laguerre(eta, n-1, 2*l+3))
    else:
        dr = 2.0 / zeta * l*eta**(l-1)* np.exp(-0.5*eta) * assoc_laguerre(eta, n, 2*l+2) + \
            eta**l * (-1.0) * np.exp(-0.5*eta) * assoc_laguerre(eta, n, 2*l+2) + \
            eta**l * np.exp(-0.5*eta) * assoc_laguerre(eta, n-1, 2*l+3)
        return prefact*dr
    
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
    try:
        return N(wigner_6j(0.5*j1,0.5*j2,0.5*j3,0.5*j4,0.5*j5,0.5*j6))
    except ValueError:
        return 0

def njs(j1, j2, j3, j4, j5, j6, j7, j8, j9):
    """
    9-j symbol
    { j1 j2 j3 }
    { j4 j5 j6 }
    { j7 j8 j9 }
    """
    try:
        return N(wigner_9j(0.5*j1, 0.5*j2, 0.5*j3, 0.5*j4, 0.5*j5, 0.5*j6, 0.5*j7, 0.5*j8, 0.5*j9))
    except ValueError:
        return 0

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
        return np.sqrt( gamma(n1+1) * gamma(n2+2*l+3) / ( gamma(n2+1) * gamma(n1+2*l+3)) ) * (4*n2+4*l+6) / (2*l+3) * 0.5*np.sqrt(2*l+1)
    if(n1==n2):
        return (4*n1+2*l+3) / (2*l+3) * 0.5 * np.sqrt(2*l+1)
    if(n1<n2):
        return np.sqrt( gamma(n2+1) * gamma(n1+2*l+3) / ( gamma(n1+1) * gamma(n2+2*l+3)) ) * (4*n1+4*l+6) / (2*l+3) * 0.5 * np.sqrt(2*l+1)

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

def Ck_one_body(j1, l1, j2,  l2, k):
    '''
    <j1||Ck||j2> = 4pi/sqrt(2k + 1) * <l1 1/2 j1 || Y^k || l1 1/2 j2>
    '''
    if (l1+l2+k)%2==1: return 0
    return thj(j1, j2, 2*k, 1, -1, 0)*np.sqrt((j1+1)*(j2+1))*((-1)**(j2/2.+k-0.5))

def YL(j1, l1, j2, l2, L):
    """
    <l1 1/2 j1 ||[Y^L] || l1 1/2 j2>
    """
    if (l1+l2+L)%2==1: return 0
    return np.sqrt((j1+1)*(j2+1)*(2*L+1))/(4*np.pi)*thj(j1, 2*L, j2, -1, 0, 1)*((-1)**(j2/2.+L-0.5))

def YLs(j1, l1, j2, l2, L, K):
    """
    <j1 ||[Y^Ls]^K || j2>
    """
    return np.sqrt(3/2)*np.sqrt((j1+1)*(j2+1)*(2*K+1))*njs(2*l1, 1, j1, l2, 1, j2, 2*L, 1, 2*K) \
        *YL(j1, l1, j2, l2, L)

def L_element(j1, l1, j2, l2):
    """
    <l1 1/2 j1 ||[L] || l1 1/2 j2>
    """
    if l1 != l2: return 0
    return np.sqrt((j1+1)*(j2+1)*l1*(l1+1)*(2*l1+1))*sjs(2*l1, 2*l1, 2, j2, j1, 1)*((-1)**(j2/2.+l1-1.5))

def Rkk_one_body(o1, o2, k):
    '''
    <j1||R(k,k)||j2> = <j1||[Y^K L]^K||j2>
    '''
    total = 0.
    for j in np.arange(np.abs(o2.j/2. - 1), o2.j/2. + 2):
        l = o2.l
        total+=Ck_one_body(o1.j, o1.l, 2*j, l, k)*sjs(2*k, 2, 2*k, o2.j, o1.j, 2*j) * \
                    np.sqrt((2*j+1)*(o2.j+1)*l*(l+1)*(2*l+1))*(-1)**(l+o2.j/2.+1.5)*sjs(2*l, 2*l, 2, o2.j, 2*j, 1)
    return total*(2*k+1)*(-1**(o1.j/2.+o2.j/2.+k))

def YL_L(j1, l1, j2, l2, L, K):
    """
    < j1 || [Y^L L]^K || j2 > 
    """
    total = 0.
    prefact = (-1)**(j1+j2+K)*np.sqrt(2*K+1)
    for j in [l2 - 1/2, l2 + 1/2]:
        total+=sjs(2*L, 2, 2*K, j2, j1, 2*j)*YL(j1, l1, 2*j, l2, L)*L_element(2*j, l2, j2, l2)
    return prefact*total 

def YL_Ls(j1, l1, j2, l2, L, K, M):
    """
    < j1 || [[Y^L L]^Ks]^M || j2 > 
    """
    prefact = (-1)**(j1+j2+M)*math.sqrt(2*M+1)
    total = 0.
    for j in [K - 1/2, K + 1/2]:
        total += sjs(2*K, 2, 2*M, j2, j1, 2*j)*YL_L(j1, l1, 2*j, l2, L, K)*math.sqrt(3/2)
    return prefact*total

def Rkk_one_body(j1, l1, j2, l2,  k):             #function overloaded for diff input formats
    '''
    <j1||R(k,k)||j2>
    '''
    total = 0.
    for j in np.arange(np.abs(j2/2. - 1), j2/2. + 1):
        total+=Ck_one_body(j1, l1, 2*j, l2, k)*sjs(2*k, 2, 2*k, j2, j1, 2*j) * \
        np.sqrt((2*j+1)*(j2+1)*l2*(l2+1)*(2*l2+1))*(-1**(l2+j2/2.+1.5))*sjs(2*l2, 2*l2, 2, j2, 2*j, 1)
    return total*(2*k+1)*(-1**(j1/2.+j2/2.+k))


def YL_L(j1, l1, j2, l2, L, K):
    """
    < j1 || [Y^L L]^K || j2 > 
    """
    total = 0.
    prefact = (-1)**((j1+j2)/2+K)*np.sqrt(2*K+1)
    for j in [l2 - 1/2, l2 + 1/2]:
        total+=sjs(2*L, 2, 2*K, j2, j1, 2*j)*YL(j1, l1, 2*j, l2, L)*L_element(2*j, l2, j2, l2)
    return prefact*total 

def YL_Ls(j1, l1, j2, l2, L, K, M):
    """
    < j1 || [[Y^L L]^Ks]^M || j2 > 
    """
    prefact = (-1)**((j1+j2)/2+M)*math.sqrt(2*M+1)
    total = 0.
    for j in [K - 1/2, K + 1/2]:
        total += sjs(2*K, 2, 2*M, j2, j1, 2*j)*YL_L(j1, l1, 2*j, l2, L, K)*math.sqrt(3/2)
    return prefact*total

def Ck_cross_s(o1, o2, k):
    '''
    <j1||[Ck x S1]k||j2>
    '''
    total = 0.
    for j in np.arange(np.abs(o2.j/2. -1), o2.j/2.+1):
        total+=sjs(2*k, 2, 2*k, o2.j, o1.j, 2*j)*Ck_one_body(o1.j, o1.l, 2*j, o2.l, k)*np.sqrt(6*(o2.j+1)*(2*j+1))*(-1**(o2.l+j+1.5))*sjs(1, 1, 2, o2.j, 2*j, 2*o2.l)
    return total*np.sqrt(2*k+1)*(-1)**(k+o1.j/2.+o2.j/2.)

def Rkk_cross_s(o1, o2, k):
    '''
    <j1j2||[Ck(1) x Ck(2)||j3j4>
    '''
    total = 0.
    for j in np.arange(np.abs(o2.j/2. -1), o2.j/2.+1):
        total+=sjs(2*k, 2, 2*k, o2.j, o1.j, 2*j)*Rkk_one_body(o1.j, o1.l, 2*j, o2.l, k)*np.sqrt(6*(o2.j+1)*(2*j+1))*(-1**(o2.l+j+1.5))*sjs(1, 1, 2, o2.j, 2*j, 2*o2.l)
    return total*np.sqrt(2*k+1)*(-1)**(k+o1.j/2.+o2.j/2.)


def Rkk1_cross_Ck2_dot_S1(o1, o2, o3, o4, J, k):
    return -Rkk_cross_s(o1, o3, k)*Ck_one_body(o2.j, o2.l, o4.j, o4.l, k)/(3*(2*k+1)*np.sqrt(2*J+1))

def Rkmin1k1_cross_Ckmin12_dot_S1(o1, o2, o3, o4, J, k):
    return -Rkk_cross_s(o1, o3, k-1)*Ck_one_body(o2.j, o2.l, o4.j, o4.l, k-1)/(3*(2*k+1)*np.sqrt(2*J+1))

def Rkpls1k1_cross_Ckpls12_dot_S1(o1, o2, o3, o4, J, k):
    return -Rkk_cross_s(o1, o3, k+1)*Ck_one_body(o2.j, o2.l, o4.j, o4.l, k+1)/(3*(2*k+1)*np.sqrt(2*J+1))


def Ck1_cross_Ck2_dot_S1(o1, o2, o3, o4, J, k):
    return -Ck_cross_s(o1, o3, k)*Ck_one_body(o2.j, o2.l, o4.j, o4.l, k)/(3*(2*k+1)*np.sqrt(2*J+1))

def Ck1_cross_Ck2_dot_S2(o1, o2, o3, o4, J, k):
    return -Ck_cross_s(o2, o4, k)*Ck_one_body(o1.j, o1.l, o3.j, o3.l, k)/(3*(2*k+1)*np.sqrt(2*J+1))

def Ck_dot_Rkk(o1, o2, o3, o4, J, k):
    return Ck_one_body(o1.j, o1.l, o3.j, o3.l, k)*Rkk_one_body(o2, o4, k)*sjs(o1.j, o2.j, 2*J, o4.j, o3.j, 2*k)*(-1**(o2.j/2.+J+o3.j/2.))


def Rkk_dot_Rkk(o1, o2, o3, o4, J, k):
    return Rkk_one_body(o1, o3, k)*Rkk_one_body(o2, o4, k)*sjs(o1.j, o2.j, 2*J, o4.j, o3.j, 2*k)*(-1**(o2.j/2.+J+o3.j/2.))

def Ck_dot_Rkmin1k(o1, o2, o3, o4, J, k):
    return Ck_one_body(o1.j, o1.l, o3.j, o3.l, k)*Rkk_one_body(o2, o4, k-1)*sjs(o1.j, o2.j, 2*J, o4.j, o3.j, 2*k)*(-1)**(o2.j/2.+J+o3.j/2.)

def Rkmin1k_dot_Rkmin1k(o1, o2, o3, o4, J, k):
    return Rkk_one_body(o1, o3, k-1)*Rkk_one_body(o2, o4, k-1)*sjs(o1.j, o2.j, 2*J, o4.j, o3.j, 2*k)*(-1)**(o2.j/2.+J+o3.j/2)

def spin_harmonic_decouple(l1, j1, l2, j2, k, L):
    sum_j = 0
    for j in np.arange(np.abs(j1-1), j1+2):
        for l in np.arange(np.abs(j-0.5), j+1.5):
            if l!=l1: continue
            if (l+L+l2)%2==1: continue
            rme = (np.sqrt(1.5*(2*j1+1)*(2*j+1))*sjs(1, 1, 2, 2*j, 2*j1,2*l)*(-1**(l+j1+1.5)))
            rme *= thj(2*j, 2*j2, 2*L, 1, -1, 0) *np.sqrt((2*L+1)*(2*j+1)*(2*j2+1))*(-1**(L-0.5+j2))
            rme *= sjs(2, 2*L, 2*k, 2*j2, 2*j1, 2*j)
            sum_j+=rme
    return sum_j * (-1**(j1+j2+k))*np.sqrt(2*L+1)


def two_body_spin_contact(o1, o2, o3, o4, J, zeta):
    lmin = max( abs(o1.j-o3.j), abs(o2.j-o4.j) )//2
    lmax = min(     o1.j+o3.j ,     o2.j+o4.j  )//2
    me = 0.
    for l in range(lmin, lmax+1):
        k_term = 0.
        for k in range(np.abs(l-1), l+2):
            k_term += sjs(o1.j, o2.j, 2*J, o4.j, o3.j, 2*k)*spin_harmonic_decouple(o1.l, o1.j/2., o3.l, o3.j/2., k, l)* \
                    spin_harmonic_decouple(o2.l, o2.j/2., o4.l, o4.j/2., k, l) * np.sqrt((2*k+1)/(3*(2*l+1)))
        me+=k_term*darwin_overlap_integral(o1.n, o1.l, o2.n, o2.l, o3.n, o3.l, o4.n, o4.l, l, zeta)[0]
    return me*(-1**(J+(o2.j+o3.j)/2))*(2/(3*c*c))


def darwin_overlap_integral(n1,l1,n2,l2,n3,l3,n4,l4,L, zeta):
    return integrate.quad(lambda r1: \
            laguerre_wave_function(r1,zeta,n1,l1) * laguerre_wave_function(r1,zeta,n2,l2) * \
            laguerre_wave_function(r1,zeta,n3,l3) * laguerre_wave_function(r1,zeta,n4,l4) * \
            r1*r1, 0, np.inf)


def two_body_darwin(o1, o2, o3, o4, J, zeta):
    r = 0.0
    lmin = max( abs(o1.j-o3.j), abs(o2.j-o4.j) )//2
    lmax = min(     o1.j+o3.j ,     o2.j+o4.j  )//2
    integral = 0.0
    for l in range(lmin,lmax+1):
        if( (o1.l+o3.l+l)%2 == 1): continue
        if( (o2.l+o4.l+l)%2 == 1): continue
        integral = darwin_overlap_integral(o1.n, o1.l, o2.n, o2.l, o3.n, o3.l, o4.n, o4.l, l, zeta)[0]
        r += integral * sjs(o1.j, o2.j, J, o4.j, o3.j, 2*l) * thj(o1.j, 2*l, o3.j, -1, 0, 1) * thj(o2.j, 2*l, o4.j, -1, 0, 1)
    return r * np.sqrt( (o1.j+1) * (o2.j+1) * (o3.j+1) * (o4.j+1) ) * (1-2*( ( (o1.j+o3.j)/2+J)%2 ) ) / (4*c*c)


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
        integral, error = coulomb_F_laguerre_wave_function(o1.n, o1.l, o2.n, o2.l, o3.n, o3.l, o4.n, o4.l, l, zeta)
        r += integral * sjs(o1.j, o2.j, 2*J, o4.j, o3.j, 2*l) * thj(o1.j, 2*l, o3.j, -1, 0, 1) * thj(o2.j, 2*l, o4.j, -1, 0, 1)
    return r * np.sqrt( (o1.j+1) * (o2.j+1) * (o3.j+1) * (o4.j+1) ) * (-1)**((o1.j+o3.j)/2+J)

def LS2_A_laguerre(n1,l1,n2,l2,n3,l3,n4,l4,L, zeta, r1, r2):
    """
    Integral A^L(a, b, c, d) in the spin-orbit coupling matrix element
    """
    if r1 >= r2:
        return laguerre_wave_function(r1, zeta, n1, l1) * laguerre_wave_function(r2, zeta, n2, l2) * \
                dlaguerre_wave_function(r1, zeta, n3, l3) * laguerre_wave_function(r2, zeta, n4, l4) * \
                r2**(L+2)/r1**(L)
    else:
        return laguerre_wave_function(r1, zeta, n1, l1) * laguerre_wave_function(r2, zeta, n2, l2) * \
                dlaguerre_wave_function(r1, zeta, n3, l3) * laguerre_wave_function(r2, zeta, n4, l4) * \
                r1**(L+1)/r2**(L-1)

def LS2_B_laguerre(n1,l1,n2,l2,n3,l3,n4,l4,L, zeta, r1, r2):
    """
    Integral B^L(a, b, c, d) in the spin-orbit coupling matrix element
    """
    radial = laguerre_wave_function(r1, zeta, n1, l1) * laguerre_wave_function(r2, zeta, n2, l2) * \
                laguerre_wave_function(r1, zeta, n3, l3) * laguerre_wave_function(r2, zeta, n4, l4)
    if r1 >= r2:
        return radial * (-L)*r2**(L+2)/r1**(L+1)
    else:
        return radial * (L+1)*r1**(L)/r2**(L-1)

def LS2_C_laguerre(n1,l1,n2,l2,n3,l3,n4,l4,L, zeta, r1, r2):
    """
    Integral C^L(a, b, c, d) in the spin-orbit coupling matrix element
    """
    radial = laguerre_wave_function(r1, zeta, n1, l1) * laguerre_wave_function(r2, zeta, n2, l2) * \
                laguerre_wave_function(r1, zeta, n3, l3) * laguerre_wave_function(r2, zeta, n4, l4)
    if r1 >= r2:
        return radial * r2**(L+2)/r1**(L+1)
    else:
        return 0

def LS2_D_laguerre(n1,l1,n2,l2,n3,l3,n4,l4,L, zeta, r1, r2):
    """
    Integral D^L(a, b, c, d) in the spin-orbit coupling matrix element
    """
    radial = laguerre_wave_function(r1, zeta, n1, l1) * laguerre_wave_function(r2, zeta, n2, l2) * \
                laguerre_wave_function(r1, zeta, n3, l3) * laguerre_wave_function(r2, zeta, n4, l4)
    if r2 >= r1:
        return radial * r1**(L+1)/r2**(L)
    else:
        return 0

def two_body_spin_orbit_(o1, o2, o3, o4, J, zeta, write=False):
    """
    Spin-orbit coupling matrix element
    write - verbose option for troubleshooting
    """
    if write:
        print('Start radial integrals')
    start_radial = time.time()

    lmin = max( abs(o1.j-o3.j), abs(o2.j-o4.j) )//2
    lmax = min(o1.j+o3.j, o2.j+o4.j  )//2
    A_L = []
    B_L = []
    C_L = []
    D_L = []

    # create and store radial integrals 
    for l in range(lmin,lmax+1):
        if( (o1.l+o3.l+l)%2 == 1) or ( (o2.l+o4.l+l)%2 == 1):
            A_L.append(0)
            B_L.append(0)
            C_L.append(0)
            D_L.append(0)
        else:
            A_L.append(integrate.dblquad(lambda r1, r2: LS2_A_laguerre(o1.n, o1.l, o2.n, o2.l, o3.n, o3.l, o4.n, o4.l, l, zeta, r1, r2), 0, np.inf, lambda r1: 0, lambda r1: np.inf)[0])
            B_L.append(integrate.dblquad(lambda r1, r2: LS2_B_laguerre(o1.n, o1.l, o2.n, o2.l, o3.n, o3.l, o4.n, o4.l, l, zeta, r1, r2), 0, np.inf, lambda r1: 0, lambda r1: np.inf)[0])
            C_L.append(integrate.dblquad(lambda r1, r2: LS2_C_laguerre(o1.n, o1.l, o2.n, o2.l, o3.n, o3.l, o4.n, o4.l, l, zeta, r1, r2), 0, np.inf, lambda r1: 0, lambda r1: np.inf)[0])
            D_L.append(integrate.dblquad(lambda r1, r2: LS2_D_laguerre(o1.n, o1.l, o2.n, o2.l, o3.n, o3.l, o4.n, o4.l, l, zeta, r1, r2), 0, np.inf, lambda r1: 0, lambda r1: np.inf)[0])
    rad_finish = time.time()
    if write:
        print(f"Radial integrals complete in {rad_finish-start_radial:.2f} seconds")
        print(f"A_L: {A_L}")
        print(f"B_L: {B_L}")
        print(f"C_L: {C_L}")
        print(f"D_L: {D_L}")
    me = 0.
    # calculate angular structure
    for L in range(lmin, lmax+1):
        coupling = sjs(o1.j, o2.j, J, o4.j, o3.j, 2*L)
        me += coupling*(math.sqrt(L*(L+1))/(2*L + 1)*A_L[L-lmin]*YLs(o1.j, o1.l, o3.j, o3.l, L, L)*YL(o2.j, o2.l, o4.j, o4.l, L) \
            + 1/(2*L + 1)*B_L[L-lmin]*YL_Ls(o1.j, o1.l, o3.j, o3.l, L, L, L)*YL(o2.j, o2.l, o4.j, o4.l, L) \
            + C_L[L-lmin]*YL(o2.j, o2.l, o4.j, o4.l, L)*((L+2)/math.sqrt((2*L+1)*(2*L+3))*YL_Ls(o1.j, o1.l, o3.j, o3.l, L, L+1,L) \
            + math.sqrt((L+1)*(L+2)/((2*L+1)*(2*L+3)))*YL_Ls(o1.j, o1.l, o3.j, o3.l, L+2, L+1, L)) \
            -0.5*math.sqrt(L*(L+1))/(2*L + 1)*A_L[L-lmin]*YL(o1.j, o1.l, o3.j, o3.l, L)*YLs(o2.j, o2.l, o4.j, o4.l, L, L) \
            +0.5/(2*L + 1)*B_L[L-lmin]*YL_L(o1.j, o1.l, o3.j, o3.l, L, L)*YLs(o2.j, o2.l, o4.j, o4.l, L, L) \
            +0.5*(L+2)*math.sqrt(1/((2*L+1)*(2*L+3)))*C_L[L-lmin]*YL_L(o1.j, o1.l, o3.j, o3.l, L, L+1)*YLs(o2.j, o2.l, o4.j, o4.l, L, L+1) \
            +0.5*math.sqrt((L+1)*(L+2)/((2*L+1)*(2*L+3)))*C_L[L-lmin]*YL_L(o1.j, o1.l, o3.j, o3.l, L+2, L+1)*YLs(o2.j, o2.l, o4.j, o4.l, L, L+1))
        
        if L > 0:
            me += coupling*(D_L[L-lmin]*math.sqrt(1/(2*L+1))*(-math.sqrt(L*(L-1)/(2*L+3))*YL_Ls(o1.j, o1.l, o3.j, o3.l, L-2, L-1, L)*YL(o2.j, o2.l, o4.j, o4.l, L) \
            -(L-1)/math.sqrt(2*L-1)*YL_Ls(o1.j, o1.l, o3.j, o3.l, L-1, L, L)*YL(o2.j, o2.l, o4.j, o4.l, L) \
            -0.5*math.sqrt(L*(L-1)/(2*L+3))*YL_L(o1.j, o1.l, o3.j, o3.l, L-2, L-1)*YLs(o2.j, o2.l, o4.j, o4.l, L, L-1) \
            -0.5*(L-1)/math.sqrt(2*L+3)*YL_L(o1.j, o1.l, o3.j, o3.l, L, L-1)*YLs(o2.j, o2.l, o4.j, o4.l, L, L-1)))
    angular_done = time.time()
    if write:
        print(f"Angular integrals complete in {angular_done-rad_finish:.2f} seconds")
    return (-1)**(o2.j/2 + o3.j/2 + J)*4*np.pi*1/(c*c)*me

def two_body_spin_orbit(o1, o2, o3, o4, J, zeta=1.0):
    return two_body_spin_orbit_(o1, o2, o3, o4, J, zeta) + \
              two_body_spin_orbit_(o2, o1, o4, o3, J, zeta)

if(__name__=="__main__"):
    #e13_a = ElectronOrbit(3, 2, 3, 23)
    #e13_b = ElectronOrbit(3, 1, 1, 16)
    #e13_c = ElectronOrbit(2, 1, 3, 12)
    #e13_d = ElectronOrbit(0, 2, 3, 8)
    # default 3 3 1 0
    e13_a = ElectronOrbit(3, 5, 11, 75)
    e13_b = ElectronOrbit(3, 5, 9, 74)
    e13_c = ElectronOrbit(1, 5, 11, 47)
    e13_d = ElectronOrbit(0, 5, 9, 35)
    J = 2
    # e13_d = ElectronOrbit(2, 0, 1, 5)
    # print(ee_laguerre_wave_function(e13_o, e13_o, e13_o, e13_o, J, 1.0))
    
    rabcd = two_body_spin_orbit(e13_a, e13_b, e13_c, e13_d, J, zeta=1.0)
    rbacd = two_body_spin_orbit(e13_b, e13_a, e13_c, e13_d, J, zeta=1.0) * (-1.0)**((e13_a.j + e13_b.j)/2 - J - 1)
    rabdc = two_body_spin_orbit(e13_a, e13_b, e13_d, e13_c, J, zeta=1.0) * (-1.0)**((e13_c.j + e13_d.j)/2 - J - 1)
    rbadc = two_body_spin_orbit(e13_b, e13_a, e13_d, e13_c, J, zeta=1.0) * (-1.0)**((e13_a.j + e13_b.j + e13_c.j + e13_d.j)/2)

    print(0.5*(rabcd + rbacd + rabdc + rbadc))
    #for x in np.arange(0.0, 10.0, 1):
    #    print(x, laguerre_wave_function(x, 1.0, 0, 0))

    #for n1 in range(5):
    #    t2 = norm_laguerre_wave_function_mom(n1,n1,l)[0]
    #    print(n1, n1, t2)
    #for l in range(5):
    #    print( laguerre_wave_function(0.0, 1.0, 3, l) )

    #print(thj(1,0,1,1,0,-1))
    #print(sjs(0,0,0,0,0,0))
    #print(coulomb_F_laguerre_wave_function(0,0, 0,0, 0,0, 0,0, 0,1))