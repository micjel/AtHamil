#!/usr/bin/env python3
def main():
    Orbits = ElectronOrbits(emax=2)
    Orbits.print_orbits()


class ElectronOrbit:
    def __init__(self, n, l, j, index):
        self.n = n
        self.l = l
        self.j = j
        self.index = index
    def get_parity(self):
        return 1 - 2*( self.l%2 )

class ElectronOrbits:
    def __init__(self,emax=0,lmax=None):
        self._nlj_to_index = {}
        self._orbits = []
        self._number_orbits = 0
        self._emax = emax
        self._lmax = emax
        if(lmax!=None): self._lmax=lmax
        index = 0
        for e in range(self._emax+1):
            for l in range( min(e,self._lmax)+1 ):
                n = e-l
                for s in [-1,1]:
                    j = 2*l+s
                    if( j < 1 ): continue
                    self._add_orbit(n,l,j,index)
                    index += 1

    def _add_orbit(self, n, l, j, index):
        self._nlj_to_index[(n,l,j)] = index
        o = ElectronOrbit(n,l,j,index)
        self._orbits.append(o)
        self._number_orbits = len(self._orbits)
    def print_orbits(self):
        print("- single-particle orbit --")
        for index in range(self.get_number_orbits() ):
            o = self.get_orbit(index)
            line = "index= {3:3d}, n={0:3d}, l={1:3d}, j={2:3d}".format(o.n, o.l, o.j, index)
            print(line)
    def get_orbits(self):
        return self._orbits
    def get_orbit(self, index):
        return self._orbits[index]
    def get_orbit_index(self, n, l, j):
        return self._nlj_to_index[(n,l,j)]
    def get_number_orbits(self):
        return self._number_orbits
    def get_emax(self):
        return self._emax
    def get_lmax(self):
        return self._lmax

if(__name__ == "__main__"):
    main()
