#!/usr/bin/env python3
from MyLibrary import *
from Orbits import ElectronOrbit, ElectronOrbits
from TwoBodySpace import TwoBodyChannel, TwoBodySpace
import numpy as np

def main():
    Orbits = ElectronOrbits(emax=0)
    TwoBody = TwoBodySpace(Orbits)
    pot = TwoBodyOperator(TwoBody)
    pot.set_operator("coulomb", 0.5)
    dar = TwoBodyOperator(TwoBody)
    dar.set_operator("darwin", 0.5)
    ssc = TwoBodyOperator(TwoBody)
    ssc.set_operator("spin_contact", 0.5)


class TwoBodyOperatorChannel:
    def __init__(self, chbra, chket):
        self.chbra = chbra
        self.chket = chket
        self.mat = np.zeros( (chbra.get_number_states(), chket.get_number_states()) )

class TwoBodyOperator:
    def __init__(self, space, rank_J=0, rank_P=1):
        self.space = space
        self.rank_J = rank_J
        self.rank_P = rank_P
        self.channels = {}
        for chbra in space.get_channels():
            for chket in space.get_channels():
                if( triangle( 2*rank_J, chbra.get_j(), chket.get_j() ) ): continue
                if( chbra.get_parity() * chket.get_parity() * rank_P == -1): continue
                self.channels[(chbra.channel_index, chket.channel_index)] = TwoBodyOperatorChannel(chbra,chket)

    def set_operator(self, opname, zeta):
        if(opname == "coulomb"):
            self._set_coulomb_term(zeta)
        if(opname == 'darwin'):
            self._set_darwin_term(zeta)
        if(opname == 'spin_contact'):
            self._set_spin_contact_term(zeta)

    def SetME(self, a, b, c, d, Jab, Jcd, ME):
        if( triangle( self.rank_J, Jab, Jcd ) ): return
        orbs = self.space.Orbits
        o1 = orbs.get_orbit(a)
        o2 = orbs.get_orbit(b)
        o3 = orbs.get_orbit(c)
        o4 = orbs.get_orbit(d)
        if( o1.get_parity() * o2.get_parity() * o3.get_parity() * o4.get_parity() * self.rank_P == -1): return
        ibra = self.space.get_channel_index( Jab, o1.get_parity()*o2.get_parity() )
        iket = self.space.get_channel_index( Jcd, o3.get_parity()*o4.get_parity() )
        chbra = self.space.get_channel( ibra )
        chket = self.space.get_channel( iket )
        bra = chbra.orbit_indices_to_index[(o1.index,o2.index)]
        ket = chket.orbit_indices_to_index[(o3.index,o4.index)]
        self.channels[(ibra,iket)].mat[bra,ket] = ME * chbra.phase[(o1.index,o2.index)] * chket.phase[(o3.index,o4.index)]

    def GetME(self, a, b, c, d, Jab, Jcd):
        if( triangle( self.rank_J, Jab, Jcd ) ): return 0.0
        orbs = self.space.Orbits
        o1 = orbs.get_orbit(a)
        o2 = orbs.get_orbit(b)
        o3 = orbs.get_orbit(c)
        o4 = orbs.get_orbit(d)
        if( o1.get_parity() * o2.get_parity() * o3.get_parity() * o4.get_parity() * self.rank_P == -1): return 0.0
        ibra = self.space.get_channel_index( Jab, o1.get_parity()*o2.get_parity() )
        iket = self.space.get_channel_index( Jcd, o3.get_parity()*o4.get_parity() )
        chbra = self.space.get_channel( ibra )
        chket = self.space.get_channel( iket )
        bra = chbra.orbit_indices_to_index[(o1.index,o2.index)]
        ket = chket.orbit_indices_to_index[(o3.index,o4.index)]
        return self.channels[(ibra,iket)].mat[bra,ket] * chbra.phase[(o1.index,o2.index)] * chket.phase[(o3.index,o4.index)]

    def _set_darwin_term(self, zeta):
        for ch in self.space.get_channels():
            j = ch.get_j()
            for bra in range(ch.get_number_states()):
                i1 = ch.index_to_orbit_index1[bra]
                i2 = ch.index_to_orbit_index2[bra]
                o1 = ch.Orbits.get_orbit(i1)
                o2 = ch.Orbits.get_orbit(i2)
                for ket in range(bra+1):
                    i3 = ch.index_to_orbit_index1[ket]
                    i4 = ch.index_to_orbit_index2[ket]
                    o3 = ch.Orbits.get_orbit(i3)
                    o4 = ch.Orbits.get_orbit(i4)

                    norm = 1.0
                    #if(i1==i2): norm /= np.sqrt(2.0)
                    if(i3==i4): norm /= np.sqrt(2.0)
                    me_abcd = two_body_darwin(o1,o2,o3,o4,j,zeta)
                    me_abdc = two_body_darwin(o1,o2,o4,o3,j,zeta) * (1-2*(( (o3.j+o4.j)/2-j-1)%2))
                    me = (me_abcd + me_abdc) * norm / np.sqrt(2.0)
                    self.SetME(i1,i2,i3,i4,j,j,me)
                    self.SetME(i3,i4,i1,i2,j,j,me)

    def _set_spin_contact_term(self, zeta):
        for ch in self.space.get_channels():
            j = ch.get_j()
            for bra in range(ch.get_number_states()):
                i1 = ch.index_to_orbit_index1[bra]
                i2 = ch.index_to_orbit_index2[bra]
                o1 = ch.Orbits.get_orbit(i1)
                o2 = ch.Orbits.get_orbit(i2)
                for ket in range(bra+1):
                    i3 = ch.index_to_orbit_index1[ket]
                    i4 = ch.index_to_orbit_index2[ket]
                    o3 = ch.Orbits.get_orbit(i3)
                    o4 = ch.Orbits.get_orbit(i4)

                    norm = 1.0
                    #if(i1==i2): norm /= np.sqrt(2.0)
                    if(i3==i4): norm /= np.sqrt(2.0)
                    me_abcd = two_body_spin_contact(o1,o2,o3,o4,j,zeta)
                    me_abdc = two_body_spin_contact(o1,o2,o4,o3,j,zeta) * (1-2*(( (o3.j+o4.j)/2-j-1)%2))
                    me = (me_abcd + me_abdc) * norm / np.sqrt(2.0)
                    self.SetME(i1,i2,i3,i4,j,j,me)
                    self.SetME(i3,i4,i1,i2,j,j,me)


    def _set_coulomb_term(self, zeta):
        for ch in self.space.get_channels():
            j = ch.get_j()
            for bra in range(ch.get_number_states()):
                i1 = ch.index_to_orbit_index1[bra]
                i2 = ch.index_to_orbit_index2[bra]
                o1 = ch.Orbits.get_orbit(i1)
                o2 = ch.Orbits.get_orbit(i2)
                for ket in range(bra+1):
                    i3 = ch.index_to_orbit_index1[ket]
                    i4 = ch.index_to_orbit_index2[ket]
                    o3 = ch.Orbits.get_orbit(i3)
                    o4 = ch.Orbits.get_orbit(i4)

                    norm = 1.0
                    #if(i1==i2): norm /= np.sqrt(2.0)
                    if(i3==i4): norm /= np.sqrt(2.0)
                    me_abcd = ee_laguerre_wave_function(o1,o2,o3,o4,j,zeta)
                    me_abdc = ee_laguerre_wave_function(o1,o2,o4,o3,j,zeta) * (1-2*(( (o3.j+o4.j)/2-j-1)%2))
                    me = (me_abcd + me_abdc) * norm / np.sqrt(2.0)
                    self.SetME(i1,i2,i3,i4,j,j,me)
                    self.SetME(i3,i4,i1,i2,j,j,me)

if(__name__ == "__main__"):
    main()
