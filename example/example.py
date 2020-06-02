#!/usr/bin/python

from sys import path
from os import getcwd

path.append("/Users/facioj/GitHub/Magnetic-TIS-natural/python/")

import slab_tools as sl

slab = sl.my_slab("out",rel=True)

slab.identify_blocks()
