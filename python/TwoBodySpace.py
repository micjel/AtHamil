#!/usr/bin/env python3
import sys
from MyLibrary import triangle
from Orbits import ElectronOrbit, ElectronOrbits
from OneBodySpace import OneBodyChannel, OneBodySpace

def main():
    Orbits = ElectronOrbits(emax=6)
    Orbits.print_orbits()
    OneBody = OneBodySpace(Orbits)
    OneBody.print_one_body_space()
    TwoBody = TwoBodySpace(Orbits)
    TwoBody.print_two_body_space()

class TwoBodyChannel:
    def __init__(self, Orbits, e2max, j, p, channel_index):
        self.Orbits = Orbits
        self.p = p
        self.j = j
        self.channel_index = channel_index
        self.index_to_orbit_index1 = []
        self.index_to_orbit_index2 = []
        self.orbit_indices_to_index = {}
        self.phase = {}
        self.number_states = 0
        index = 0
        for o1 in self.Orbits.get_orbits():
            for o2 in self.Orbits.get_orbits():
                if( o2.index > o1.index ): continue
                if( o1.index==o2.index and j%2==1): continue
                if( 1-2*((o1.l+o2.l)%2) != p ): continue
                if( triangle(o1.j, o2.j, 2*j) ): continue
                if( o1.n+o1.l + o2.n+o2.l > e2max): continue
                self.index_to_orbit_index1.append(o1.index)
                self.index_to_orbit_index2.append(o2.index)
                self.orbit_indices_to_index[(o1.index, o2.index)] = index
                self.orbit_indices_to_index[(o2.index, o1.index)] = index
                self.phase[(o1.index,o2.index)] = 1
                self.phase[(o1.index,o2.index)] = 1-2*( ((o1.j+o2.j)/2 - j)%2 )
                index += 1
        self._number_states = len(self.index_to_orbit_index1)
    def get_j(self):
        return self.j
    def get_parity(self):
        return self.p
    def get_number_states(self):
        return self._number_states

class TwoBodySpace:
    def __init__(self, Orbits, e2max=None):
        self.Orbits = Orbits
        self._number_channels = 0
        self._channels = []
        self._jp_to_index = {}
        self.jmax = 2*Orbits.get_lmax()+1
        self.e2max = e2max
        if(e2max == None): self.e2max=2*Orbits.get_emax()
        index = 0
        for j in range(self.jmax+1):
            for p in [1,-1]:
                chan = TwoBodyChannel(Orbits, self.get_e2max(), j, p, index)
                if( chan.get_number_states() == 0): continue
                self._jp_to_index[(j,p)] = index
                self._channels.append(chan)
                index += 1
        self._number_channels = len(self._channels)
    def get_channel(self,index):
        return self._channels[index]
    def get_channels(self):
        return self._channels
    def get_channel_index(self,j,p):
        return self._jp_to_index[(j,p)]
    def get_number_channels(self):
        return self._number_channels
    def get_emax(self):
        return self.Orbits.get_emax()
    def get_e2max(self):
        return self.e2max
    def get_lmax(self):
        return self.Orbits.get_lmax()
    def print_two_body_space(self):
        print("- Two-body space --")
        for ch in range(self.get_number_channels()):
            chan = self.get_channel(ch)
            line = "index= {0:3d}, j={1:3d}, p={2:3d}, # of states={3:6d}".format(ch, chan.get_j(), chan.get_parity(), chan.get_number_states())
            print(line)



if(__name__ == "__main__"):
    main()
