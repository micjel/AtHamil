#!/usr/bin/env python3
import sys
from Orbits import ElectronOrbit, ElectronOrbits
def main():
    Orbits = ElectronOrbits(emax=6)
    Orbits.print_orbits()
    OneBody = OneBodySpace(Orbits)
    OneBody.print_one_body_space()

class OneBodyChannel:
    def __init__(self, Orbits, j, p, channel_index):
        self.Orbits = Orbits
        self.p = p
        self.j = j
        self.channel_index = channel_index
        self.index_to_orbit_index = []
        self.orbit_index_to_index = {}
        self.number_states = 0
        index = 0
        for o in self.Orbits.get_orbits():
            if(1-2*(o.l%2) != p): continue
            if(o.j != j): continue
            self.index_to_orbit_index.append(o.index)
            self.orbit_index_to_index[o.index] = index
            index += 1
        self._number_states = len(self.index_to_orbit_index)
    def get_j(self):
        return self.j
    def get_parity(self):
        return self.p
    def get_number_states(self):
        return self._number_states

class OneBodySpace:
    def __init__(self, Orbits):
        self.Orbits = Orbits
        self._number_channels = 0
        self._channels = []
        self._jp_to_index = {}
        self.jmax = 2*Orbits.get_lmax()+1
        index = 0
        for j in range(1, self.jmax+2, 2):
            for p in [1,-1]:
                chan = OneBodyChannel(Orbits, j, p, index)
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
    def get_lmax(self):
        return self.Orbits.get_lmax()
    def print_one_body_space(self):
        print("- One-body space --")
        for ch in range(self.get_number_channels()):
            chan = self.get_channel(ch)
            line = "index= {0:3d}, j={1:3d}, p={2:3d}, # of states={3:3d}".format(ch, chan.get_j(), chan.get_parity(), chan.get_number_states())
            print(line)

if(__name__ == "__main__"):
    main()
