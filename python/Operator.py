#!/usr/bin/env python3
from MyLibrary import triangle, T_laguerre_wave_function, Coul_laguerre_wave_function1
from ModelSpace import ModelSpace
from OneBodyOperator import OneBodyOperator
from TwoBodyOperator import TwoBodyOperator
import numpy as np

def main():
    ms = ModelSpace(emax=1, zeta=2.0)
    Ham = Operator(ms)
    #Ham.set_operator("kinetic")
    Ham.set_operator("coulomb")
    for chan in Ham.OneBody.channels.values():
        print(chan.mat)
    for chan in Ham.TwoBody.channels.values():
        print(chan.mat)

class Operator:
    def __init__(self, ms):
        self.ms = ms
        self.OneBody = OneBodyOperator( ms.get_one_body() )
        self.TwoBody = TwoBodyOperator( ms.get_two_body() )
    def set_operator(self, opname):
        zeta_inv = 1.0/self.ms.zeta
        self.OneBody.set_operator(opname, zeta_inv)
        self.TwoBody.set_operator(opname, zeta_inv)


if(__name__ == "__main__"):
    main()
