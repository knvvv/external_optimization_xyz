from relaxpack.graphbuilder import IK_Molecule,drawgraph
from ctypes import *
import os
mol = IK_Molecule("mysubj","molfromsmiles.sdf")
mol.createcpp()
mol.compile()
N = len(mol.atoms_xyz)

class molgeometry(Structure):
    _fields_ = [('nconf', c_int), ('arr', c_double*(3*N)), ('symbols', c_char_p*(N))]


print(repr(os.environ['LD_LIBRARY_PATH']))

so_file = "./%s.so" % "mylib" #mol.name
dll = CDLL(so_file)
funct = dll.relaxgeom
funct.argtypes = [POINTER(molgeometry)]
# RANDOM COMMENT
symbols = []
mg = molgeometry()
xyz = open("source.xyz", "r").readlines()
for i in range(N):
    parts = list(filter(None, xyz[i].split(" ")))
    symbols.append(parts[0])
    mg.symbols[i] = str.encode(parts[0][0])

    mg.arr[3*i] = float(parts[1])
    mg.arr[3*i+1] = float(parts[2])
    mg.arr[3*i+2] = float(parts[3])
funct(byref(mg))
print("I'm alive")

geom = []
for i in range(N):
    geom.append([mg.arr[3*i],mg.arr[3*i+1],mg.arr[3*i+2]])
print("I'm alive2")
xyz = open("relaxedgeom.xyz","w")
for i in range(N):
    xyz.write("%s %f %f %f\n" % (symbols[i],geom[i][0],geom[i][1],geom[i][2]))
xyz.close()
print("I'm alive3")


