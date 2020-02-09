#!/usr/bin/env python3
import sys
from Orbits import ElectronOrbit, ElectronOrbits
from OneBodySpace import OneBodyChannel, OneBodySpace
from TwoBodySpace import TwoBodyChannel, TwoBodySpace

def main():
    ms = ModelSpace()

class ModelSpace:
    def __init__(self, emax=6, e2max=None, lmax=None, zeta=1.0):
        self.zeta = zeta
        self.emax = emax
        self.e2max = e2max
        self.lmax = lmax
        if(e2max==None): self.e2max = 2*self.emax
        if(lmax==None): self.lmax = 2*self.emax
        self.Orbits = ElectronOrbits(self.emax, lmax=self.lmax)
        self.OneBody = OneBodySpace(self.Orbits)
        self.TwoBody = TwoBodySpace(self.Orbits, e2max=self.e2max)
    def get_emax(self):
        return self.emax
    def get_e2max(self):
        return self.e2max
    def get_lmax(self):
        return self.lmax
    def get_one_body(self):
        return self.OneBody
    def get_two_body(self):
        return self.TwoBody

if(__name__ == "__main__"):
    main()
