#!/usr/bin/env python3
from MyLibrary import triangle, T_laguerre_wave_function, Coul_laguerre_wave_function1
from Orbits import ElectronOrbit, ElectronOrbits
from OneBodySpace import OneBodyChannel, OneBodySpace
import numpy as np

def main():
    Orbits = ElectronOrbits(emax=2)
    OneBody = OneBodySpace(Orbits)
    kin = OneBodyOperator(OneBody)
    kin.set_operator("kinetic",2)
    pot = OneBodyOperator(OneBody)
    pot.set_operator("coulomb",2)

class OneBodyOperatorChannel:
    def __init__(self, chbra, chket):
        self.chbra = chbra
        self.chket = chket
        self.mat = np.zeros( (chbra.get_number_states(), chket.get_number_states()) )

class OneBodyOperator:
    def __init__(self, space, rank_J=0, rank_P=1):
        self.space = space
        self.rank_J = rank_J
        self.rank_P = rank_P
        self.channels = {}
        for chbra in space.get_channels():
            for chket in space.get_channels():
                if( triangle( 2*rank_J, chbra.get_j(), chket.get_j() ) ): continue
                if( chbra.get_parity() * chket.get_parity() * rank_P == -1): continue
                self.channels[(chbra.channel_index, chket.channel_index)] = OneBodyOperatorChannel(chbra,chket)

    def set_operator(self, opname, zeta):
        if(opname == "kinetic"):
            self._set_kinetic_term(zeta)
        if(opname == "coulomb"):
            self._set_coulomb_term(zeta)
    def _set_kinetic_term(self, zeta):
        orbs = self.space.Orbits
        for o1 in orbs.get_orbits():
            for o2 in orbs.get_orbits():
                if( o1.j != o2.j  ): continue
                if( o1.l != o2.l  ): continue
                me = T_laguerre_wave_function(o1.n, o2.n, o1.l) / zeta**2
                self.SetME( o1.index, o2.index, me )
    def _set_coulomb_term(self, zeta):
        orbs = self.space.Orbits
        for o1 in orbs.get_orbits():
            for o2 in orbs.get_orbits():
                if( o1.j != o2.j  ): continue
                if( o1.l != o2.l  ): continue
                me = Coul_laguerre_wave_function1(o1.n, o2.n, o1.l,zeta)
                self.SetME( o1.index, o2.index, me[0] )

    def SetME(self, a, b, ME):
        orbs = self.space.Orbits
        o1 = orbs.get_orbit(a)
        o2 = orbs.get_orbit(b)
        if( triangle( 2*self.rank_J, o1.j, o2.j ) ): return
        if( o1.get_parity() * o2.get_parity() * self.rank_P == -1): return
        ibra = self.space.get_channel_index( o1.j, o1.get_parity() )
        iket = self.space.get_channel_index( o2.j, o2.get_parity() )
        chbra = self.space.get_channel( ibra )
        chket = self.space.get_channel( iket )
        bra = chbra.orbit_index_to_index[o1.index]
        ket = chket.orbit_index_to_index[o2.index]
        self.channels[(ibra,iket)].mat[bra,ket] = ME

    def GetME(self, a, b):
        orbs = self.space.Orbits
        o1 = orbs.get_orbit(a)
        o2 = orbs.get_orbit(b)
        if( triangle( 2*self.rank_J, o1.j, o2.j ) ): return 0
        if( o1.get_parity() * o2.get_parity() * self.rank_P == -1): return 0
        ibra = self.space.get_channel_index( o1.j, o1.get_parity() )
        iket = self.space.get_channel_index( o2.j, o2.get_parity() )
        chbra = self.space.get_channel( ibra )
        chket = self.space.get_channel( iket )
        bra = chbra.orbit_index_to_index[o1.index]
        ket = chket.orbit_index_to_index[o2.index]
        return self.channels[(ibra,iket)].mat[bra,ket]

if(__name__ == "__main__"):
    main()
