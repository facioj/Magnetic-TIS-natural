#!/usr/bin/python

import sys
file = open(sys.argv[1])
data = file.readlines()
file.close()

fila0 = data[0]
vals = map(eval,fila0[1:].split())
Ne = vals[2]
Np = vals[4]
print Ne,Np
file = open("bands_full.dat",'w')


for I in range(Np):
  for J in range(Ne):
      line = I + Np * J
      fila = data[2+line]
      vals = map(eval,fila.split())
      #print >> file, vals[0],vals[1],vals[2],vals[3],vals[4],vals[5]
      for col in range(len(vals)):
          print >> file, vals[col],
      print >> file,""
          #vals[1],vals[2],vals[3],vals[4],vals[5]
  print >> file, "\n",
file.close()
