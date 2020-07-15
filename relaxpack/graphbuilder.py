import networkx as nx
from copy import deepcopy,copy
import time
import pylab
import subprocess
import matplotlib.pyplot as plt
import sympy as sp
import numpy as np
from numpy.linalg import norm
from sympy.core.evaluate import evaluate

def drawgraph(G):
    edge_colors = ['black' for edge in G.edges()]
    pos=nx.spring_layout(G)
    nx.draw(G,pos, with_labels=True,node_size=200,edge_color=edge_colors,edge_cmap=plt.cm.Reds)
    pylab.show()

def runc(runline):
    p = subprocess.Popen(runline, shell = True)
    while p.poll() == None:
        time.sleep(2)

class IK_Molecule:
    def __init__(self, name, sdffile):
        self.name = name
        self.atoms_xyz = []
        self.atoms_sym = []
        self.bonds_num = []
        self.bonds_type = []
        self.readsdf(sdffile)
        self.creategraph()

        self.x = []
        self.y = []
        self.z = []
        self.createsymbols()
        self.geterf()

        self.grad = []
        self.getgrad()

    def readsdf(self, file):
        lines = open(file, "r").readlines()
        headdata = lines[3].replace("\n", "").split(" ")
        headdata = list(filter(None, headdata))
        try:
            natoms = int(headdata[0])
            nbonds = int(headdata[1])
        except:
            raise Exception("Invalid molfile format")
        for i in range(4, 4 + natoms):
            parts = lines[i].replace("\n", "").split(" ")
            parts = list(filter(None, parts))
            self.atoms_xyz.append(np.array([float(parts[0]), float(parts[1]), float(parts[2])]))
            self.atoms_sym.append(parts[3])
        for i in range(4 + natoms, 4 + natoms + nbonds):
            parts = lines[i].replace("\n", "").split(" ")
            parts = list(filter(None, parts))
            self.bonds_num.append([int(parts[0]) - 1, int(parts[1]) - 1])
            self.bonds_type.append(int(parts[2]))
        print("File is read")

    def creategraph(self):
        self.G = nx.Graph()
        for bond in self.bonds_num:
            self.G.add_edge(bond[0],bond[1])

    def createsymbols(self):
        N = len(self.atoms_xyz)
        for i in range(N):
            self.x.append(sp.Symbol("coord["+str(3*i)+"]"))
            self.y.append(sp.Symbol("coord["+str(3*i+1)+"]"))
            self.z.append(sp.Symbol("coord["+str(3*i+2)+"]"))

    def getpairs(self,node):
        res = []
        neighb = list(self.G.neighbors(node))
        for i in range(len(neighb)):
            for j in range(i+1,len(neighb)):
                res.append([neighb[i],neighb[j]])
        return res

    def getdot(self,node,pair):
        return np.dot(self.atoms_xyz[pair[0]]-self.atoms_xyz[node], self.atoms_xyz[pair[1]]-self.atoms_xyz[node])

    def geterf(self):
        self.erf = 0
        with evaluate(False):
            for edge in list(self.G.edges):
                self.erf += ((self.x[edge[1]] - self.x[edge[0]])**2 + (self.y[edge[1]] - self.y[edge[0]])**2 +
                        (self.z[edge[1]] - self.z[edge[0]])**2 - norm(self.atoms_xyz[edge[0]]-self.atoms_xyz[edge[1]]))**2

        for node in list(self.G.nodes):
            for pair in self.getpairs(node):
                self.erf += ( (self.x[pair[0]]-self.x[node])*(self.x[pair[1]]-self.x[node]) + (self.y[pair[0]]-self.y[node])*(self.y[pair[1]]-self.y[node]) +
                         (self.z[pair[0]]-self.z[node])*(self.z[pair[1]]-self.z[node]) - self.getdot(node,pair) )**2

    def getgrad(self):
        for i in range(len(self.atoms_xyz)):
            self.grad.append(sp.diff(self.erf, self.x[i]))
            self.grad.append(sp.diff(self.erf, self.y[i]))
            self.grad.append(sp.diff(self.erf, self.z[i]))
            print("Done with %d atom" % i)

        for item in self.grad:
            print(sp.ccode(item))

    def getdistatoms(self):
        pairs = []
        nster=0
        for i in range(len(self.grad)):
            for j in range(i,len(self.atoms_xyz)):
                if nx.shortest_path_length(self.G,source=i,target=j) > 2:
                    pairs.append("distatoms[0][%d] = %d;distatoms[1][%d] = %d;" % (nster, i, nster, j))
                    nster += 1
        self.NStericPairs = len(pairs)
        return "".join(pairs)

    def createcpp(self):
        errline = sp.ccode(self.erf)
        #apply("return %s;" % (errline), "$ERRF$")

        gradlines = []
        for i in range(len(self.grad)):
            gradlines.append("mygrad[%d] = %s;" % (i, sp.ccode(self.grad[i])))

        cpplines = open("cfiles/cycledef_template.cpp", "r").readlines()
        for i in range(len(cpplines)):
            if "$GRAD$" in cpplines[i]:
                cpplines[i] = cpplines[i].replace("$GRAD$","\n".join(gradlines))
            elif "$ERRF$" in cpplines[i]:
                cpplines[i] = cpplines[i].replace("$ERRF$", "return %s+getStericEnergy();" % errline)
            elif "$ATCONST$" in cpplines[i]:
                cpplines[i] = cpplines[i].replace("$ATCONST$", self.getdistatoms())

        newcpp = open("cfiles/cycledef.cpp","w")
        newcpp.write("".join(cpplines))
        newcpp.close()

        hlines = open("cfiles/mainrand_template.h", "r").readlines()
        for i in range(len(hlines)):
            if "$NUMATS$" in hlines[i]:
                hlines[i] = hlines[i].replace("$NUMATS$", str(len(self.atoms_xyz)))
            elif "$ATNUM$" in hlines[i]:
                hlines[i] = hlines[i].replace("$ATNUM$", str(self.NStericPairs))
        newh = open("cfiles/mainrand.h", "w")
        newh.write("".join(hlines))
        newh.close()

    def compile(self):
        runc("make")