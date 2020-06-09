#!/usr/bin/python

from sys import path
from os import getcwd

path.append("/home/facioj/Magnetic-TIS-natural/python/")

import slab_tools as sl

slab = sl.my_slab("out",rel=True)

slab.identify_blocks()
#slab.orbital_exp_decay(lambda_0=10./sl.to_angs,only_Bi=False,only_Te=False)


sites = [41,17,45,15,39,31,7,23,1,21,5,29]
orbitals = ['6p3/2-3/2','6p3/2+3/2','6p3/2-1/2','6p3/2+1/2','6p1/2-1/2','6p1/2+1/2']
orbitals = ['5p3/2-3/2','5p3/2+3/2','5p3/2-1/2','5p3/2+1/2','5p1/2-1/2','5p1/2+1/2']
orbitals = ['3d+1up']
sites = [1]
slab.orbital_dos_exp_decay(sites,orbitals,lambda_0=10./sl.to_angs,dos_folder="./dos",type_file="l")

