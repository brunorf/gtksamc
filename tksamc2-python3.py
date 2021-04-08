#!/usr/bin/env python
# coding: utf-8

__description__ = \
"""
GTKSA - Electrostatic Free Energy calculation for each ionizable residue
"""

__author__ = "Vinicius G. Contessoto, Guilherme V Bossa, Vinicius M Oliveira"
__date__ = "21/12/2016"

################################################################
#
# Version 2.0
#
# python tksamc.py -h # for help
# 
# The following programs are provided free of charge only for academic use. 
# By downloading these programs you implicitly agree that they will be used exclusively in academic research.
#
################################################################

import sys
import numpy as np
import scipy as sc
import pandas as pd
import math3d as m3d
import sympy
import subprocess
import os
import argparse
import time
from mpl_toolkits.mplot3d import Axes3D
from itertools import chain
from scipy.spatial import distance
from scipy.special import eval_legendre
from numpy import linalg as LA
from scipy.spatial.transform import Rotation as R
from itertools import islice
from subprocess import call
import profile
import threading
import nestle
import mdtraj as md

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pylab import *



parser = argparse.ArgumentParser(description='Charge-charge energy calculation in python')
parser.add_argument('-ph', action='store', default=7.0, dest='arg_pH', help='pH value')              
parser.add_argument('-T', action='store', default=300.0, dest='arg_T',  help='Temperature value')           
parser.add_argument('-f', metavar='input-file-PDB',help='insert a PDB file',type=argparse.FileType('rt'))
parser.add_argument('-e', action='store',choices=['TK'], default="TK",dest='arg_e',type=str,help='Electrostatic energy calculation method')
parser.add_argument('-s', action='store',choices=['EX','MC'], default="MC",dest='arg_s',type=str,help='Statistical method to protonation state amostration - EX = Exact; MC = Monte Carlo;')
parser.add_argument('-plot', action='store',choices=['yes','no'], default="yes",dest='arg_plot',type=str,help='Save Plot figure file - EPS')



try:
    arguments = parser.parse_args()
    print('################################################')
    print(u"\U0001F63A", "### TKSA started ###", u"\U0001F63A")
    print('Input file:', arguments.f.name)
    print('pH  =', arguments.arg_pH)
    print('T   =', arguments.arg_T)
    print('Elec. Energy Calc Method =', arguments.arg_e)
    print('Statistical Method =', arguments.arg_s)
    print('Plot =', arguments.arg_plot)
    
except IOError:
    parser.error() 


All_residues = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
Area_residues = [113,140,151,183,218,85,194,182,211,180,204,158,143,189,241,122,146,160,259,229] # colocar a referencia
Charged_residues = ['ARG','LYS','N_TER','HIS','GLU','ASP','C_TER']
Charge_values = [0,0,0,0,-1,-1,-1]
Charged_atoms = ['NH2','NZ','NE2','OE2','OD2']
PKA = [12.0,10.6,7.7,6.3,4.5,4.0,3.6]
e = (4.0,78.5,1.0,6.02,1.602,8.8541878176,1.3806488,8.314) #(ep,es,x,mol,e,eo,K,R) reduzidos



def rot_matrix(X,Y,Z):
## This function rotate all coordinates and let the atoms that present the largest distance between them in Z = 0 ##
    X = X - np.mean(X)
    Y = Y - np.mean(Y)
    Z = Z - np.mean(Z)
    
    XYZ = np.vstack((X,Y,Z)).T
    Origin = np.zeros(np.shape(XYZ))
    dist = distance.cdist(XYZ,XYZ, 'euclidean')
    resid_max_dist = np.where(dist == dist.max())[0]
    resid_a = resid_max_dist[0]
    resid_b = resid_max_dist[1]
    max_dist = np.max(dist)   

    X = X - X[resid_a] 
    Y = Y - Y[resid_a]
    Z = Z - Z[resid_a]
  

    #calculate rotation angles

    thetaY = np.arctan( X[resid_b]/Y[resid_b])

    Xn = X*np.cos(thetaY) - Y*np.sin(thetaY)
    Yn = X*np.sin(thetaY) + Y*np.cos(thetaY)

    thetaZ = np.arctan(Yn[resid_b]/Z[resid_b])
    Ynn = Yn*np.cos(thetaZ) - Z*np.sin(thetaZ)
    Zn = Yn*np.sin(thetaZ) + Z*np.cos(thetaZ)

    ## Y must be larger than X
    aa = []
    bb = []
    for i in range(len(Ynn)):
       for j in range(len(Ynn)):
          aa = np.append(aa, (np.subtract(Ynn[i], Ynn[j])))
          bb = np.append(bb, (np.subtract(Xn[i], Xn[j])))
    
    difY = np.amax(aa)
    difX = np.amax(bb)

    
    ## X must be larger than Y and X>Y>Z
    ## Invert coordinates X and Z

    Znn = Xn
    Xnn = Zn
    
    if difX > difY:
        Znn = Ynn
        Ynn = Xn

#    back to CM in X and Y coordinates
    Xnn = Xnn - np.mean(Xnn)
    Ynn = Ynn - np.mean(Ynn)
    Znn = Znn - np.mean(Zn)
#    Znn = Zn - (Zn[resid_b]/2.0)
    XYZn = np.vstack((Xnn,Ynn,Znn)).T
    return(XYZn)

##################################################################################################
# # Lame function of first kind
##################################################################################################
def lameE(aaa,bbb,ccc):
    if (aaa == 0 and bbb == 1):
        return(1.0)
    if aaa == 1 and bbb == 1:
        return(ccc)
    if aaa == 1 and bbb == 2:
        return(np.power((ccc**2 - h**2),0.5))
    if aaa == 1 and bbb == 3:
        return(np.power((ccc**2 - k**2),0.5))
                
##################################################################################################
# # Lame function of second kind
##################################################################################################
def lameF(nn,pp,ccc):
    return(lameE(nn,pp,ccc)/(ccc**(2*nn+1)))


##################################################################################################
# # Kirkwood Polynomial Function
##################################################################################################
def Kn(n,x):
   
   Kpf = np.sum([np.power(x,s)*np.divide(np.power(2.0,s)*np.math.factorial(n)*np.math.factorial(2*n-s),np.math.factorial(s)*np.math.factorial(2*n)*np.math.factorial(n-s)) for s in range(n)])
   return Kpf
##################################################################################################


#global Q,E,S,Pk,e,T,pH,total_charged_residues,G,G0,indiv_data,Gqq,Q0

file_pdb = arguments.f     
pH = np.float(arguments.arg_pH)
T = np.float(arguments.arg_T)


##################################################################################################
   # runs the standalone version of ©Surfrace
##################################################################################################

      
print('Running SASA - ©Surfrace')
cmd1 = 'echo 1' + arguments.f.name + ' 1.4 1| ./surfrace5_0_linux_64bit > SASA_'+os.path.splitext(arguments.f.name)[0]+'_all.trash' ## Roda o programa para a SASA
os.system(cmd1)
try: 
  file_sasa = open(os.path.splitext(arguments.f.name)[0] + '_residue.txt', 'r') ## Abre o arquivo que vem do programa acima
except (IOError) as errno:
  print ('I/O error - ** Check the files of SASA calculation - something went wrong **. %s' % errno)
  sys.exit()
   
SASA_data=[]
for line2 in file_sasa:
   list2 = line2.split()
   Area_norm = np.float(list2[2])/np.float(Area_residues[All_residues.index(list2[1])])
   if Area_norm >= 1.0:
       print("Warning - ** SASA greater than 1.0**",list2[1],list2[0],list2[2],np.float(Area_residues[All_residues.index(list2[1])]),Area_norm)
       print("Automatically changed to 0.75")
       Area_norm = 0.750000000001
   SASA_data.append([list2[1],list2[2],Area_norm])
indiv_data=[]
S=[]
SAij=[]
total_atoms=[]
total_residues=[]
total_charged_residues=[]


for line in file_pdb: ## Reading file.pdb
 lista = line.split()
 id = lista[0]
 if id == 'ATOM':
   atom_index = np.int(lista[1]) 
   atom_type = lista[2]
   residue_type = lista[3]
   chain = lista[4]
   residue_index = np.int(lista[5])
   total_atoms.append([atom_index])
   if atom_type == 'CA' and chain == 'A':
      total_residues.append([residue_index])
   if atom_index == 1 and atom_type == 'N' and chain == 'A' and residue_index == 1 and not residue_type in Charged_residues: ## Select the charged residues
      total_charged_residues.append([atom_index])
      S.append(['N_T',residue_index,np.size(total_charged_residues),lista[1],lista[2],np.float(lista[6]),np.float(lista[7]),np.float(lista[8]),PKA[Charged_residues.index('N_TER')],SASA_data[np.size(total_residues)-1][2],Charge_values[Charged_residues.index('N_TER')]])
   if residue_type in Charged_residues and atom_type in Charged_atoms: ## Seleciona os resíduos carregados
      total_charged_residues.append([atom_index])
      S.append([lista[3],residue_index,np.size(total_charged_residues),lista[1],lista[2],np.float(lista[6]),np.float(lista[7]),np.float(lista[8]),PKA[Charged_residues.index(residue_type)],SASA_data[np.size(total_residues)-1][2],Charge_values[Charged_residues.index(residue_type)]])
   if atom_type == 'OXT' and chain == 'A'  and not residue_type in Charged_residues:
      total_charged_residues.append([atom_index])
      S.append(['C_T',residue_index,np.size(total_charged_residues),lista[1],lista[2],np.float(lista[6]),np.float(lista[7]),np.float(lista[8]),PKA[Charged_residues.index('C_TER')],SASA_data[np.size(total_residues)-1][2],Charge_values[Charged_residues.index('C_TER')]])

print("There are: %d Charged_residues" % np.size(total_charged_residues))
Restype=np.asarray([i[0] for i in S])
X=np.asarray([i[5] for i in S])
Y=np.asarray([i[6] for i in S])
Z=np.asarray([i[7] for i in S])
XYZold = np.vstack((X,Y,Z)).T

np.savetxt('x.txt',X)   
np.savetxt('y.txt',Y) 
np.savetxt('z.txt',Z)

XYZ = rot_matrix(X,Y,Z) # Rotate coordinates for Z>Y>X
dist = distance.cdist(XYZ,XYZ, 'euclidean')

X=np.asarray([i[0] for i in XYZ])
Y=np.asarray([i[1] for i in XYZ])
Z=np.asarray([i[2] for i in XYZ])

resid_max_dist = np.where(dist == dist.max())[0]
resid_a = resid_max_dist[0]
resid_b = resid_max_dist[1]
max_dist = np.max(dist)   
Origin = np.zeros(np.shape(XYZ))
#theta_origin = np.arccos(1-angle_origin)
dist_origin = distance.cdist(XYZ, Origin, 'euclidean')
angle = distance.cdist(XYZ, XYZ, 'cosine')
raio = (np.max(dist)*0.5 + 3.4+2.0, np.max(dist)*0.5 + 2.0+2.0)
np.seterr(invalid='ignore')
np.seterr(divide='ignore')
#theta = np.arccos(1-angle)
NormA = np.matrix([LA.norm(v) for v in np.array(XYZ)])
rirj = np.array(np.dot(np.transpose(NormA),NormA))

tabledist = []
maxX = []
maxY = []
maxZ = []
# Max distances in each axis 
#Zmax = 2*np.max(dist_origin)
for i in range(len(Y)):
    for j in range(len(Y)):
#    tabledist = np.append(tabledist,[dist_origin[i][i], i, np.int(len(Y))])
       maxX = np.append(maxX, (np.subtract(X[i], X[j])))
       maxY = np.append(maxY, (np.subtract(Y[i], Y[j])))
       maxZ = np.append(maxZ, (np.subtract(Z[i], Z[j])))

Xmax = np.amax(maxX)
Ymax = np.amax(maxY)
Zmax = np.amax(maxZ)

#tabledist = tabledist.reshape((len(Y), 3))
lc = [Xmax,Ymax,Zmax]
ver=1

#optimizing axis
cr = 20 # initial high number for eccentricity 
tcheck = [0]
while (np.max(tcheck)<=0.99):
    r1=(Xmax*0.5)+(cr)
    r2=(Ymax*0.5)+(cr)
    r3=(Zmax*0.5)+(cr)
    tcheck=[(X/r1)**2+(Y/r2)**2+(Z/r3)**2];
    cr = cr - 0.01
    
# Determining semi-axis of the ellipsoid
aa=r1 + 2
bb=r2 + 2
cc=r3 + 2
h=np.power((aa**2 - bb**2),0.5)
k=np.power((aa**2 - cc**2),0.5)

w1 = -(X**2 + Y**2 + Z**2 + h**2 + k**2)
w2 = (X**2) * (h**2 + k**2) + (Y**2 * k**2) + (Z**2 * h**2) + (h**2 * k**2)
w3 = -(X**2 * h**2 * k**2)
Q = (w1**2 - 3*w2)/9
R = ((9*w1 * w2) - (27 * w3) - (2 * w1**3))/54
theta = np.arccos(R / np.power((Q**3),0.5))

s1 = np.sign(X) * np.sign(Y) * np.sign(Z)
s2 = np.sign(X) * np.sign(Y)
s3 = np.sign(X) * np.sign(Z)

# These are coordinates for the ellipsoid
lambda1 = s1 * np.power((2*np.power(Q,0.5) * np.cos(theta/3) - w1/3),0.5)
mu = s2 * np.power((2*np.power(Q,0.5) * np.cos((theta/3) + (4 * np.pi/3)) - w1/3),0.5)
nu = s3 * np.power(abs((2*np.power(Q,0.5) * np.cos((theta/3) + (2 * np.pi/3)) - w1/3)),0.5)

Xell = np.power((lambda1**2 * mu**2 * nu**2 / (k**2 * h**2)),0.5)
Yell = np.power(((lambda1**2 - h**2) * (mu**2 - h**2) * (h**2 - nu**2) / (h**2 * (k - h)**2)),0.5)
Zell = np.power(((lambda1**2 - k**2) * (k**2 - mu**2) * (k**2 - nu**2) / (k**2 * (k - h)**2)),0.5)


# Lame functions defined

#(*-with epsi the dielectric constant inside and epso outside the ellipsoid--*)
gamma=(e[0]-e[1])/(e[0]+e[1])

qind  = gamma * lameF(0,1,lambda1/lameF(0,1,aa))

xind = (lameF(0,1,aa)/lameF(0,1,lambda1)) * (lameF(1,1,lambda1 * lameE(1,1,aa))) / (lameF(1,1,aa) * lameE(1,1,lambda1)) * X
yind = (lameF(0,1,aa)/lameF(0,1,lambda1)) * (lameF(1,2,lambda1 * lameE(1,2,aa))) / (lameF(1,2,aa) * lameE(1,2,lambda1)) * Y
zind = (lameF(0,1,aa)/lameF(0,1,lambda1)) * (lameF(1,3,lambda1 * lameE(1,3,aa))) / (lameF(1,3,aa) * lameE(1,3,lambda1)) * Z
XYZind = np.vstack((xind,yind,zind)).T
distInd = distance.cdist(XYZ,XYZind, 'euclidean')


ra = r1
#rb = r1 + 2.0
#(*----The second-order image approxiamtion introduced by Deng &Cai, Comm. Comp. Phys., 2007----*)
#(*-------ka is the inverse Debye length; it is equal to 1nm for 100mM of 1:1 electrolyte-------*)
u = e[2] * aa
sig = (1-gamma) / 2.0
del1 = gamma * (1+gamma) / 2.0

#(*--Cartesian distance between target charge j and source charge i---*)
potRF = []
potC = []
pBorn = []
pCor2 = []
for i in range(len(Y)):
    for j in range(len(Y)):
        potRF = np.append(potRF, qind[i]/(4* np.pi * e[0] * distInd[i][j]))
        potC = np.append(potC, 1/(4* np.pi * e[0] * dist[i][j]))
        pBorn = np.append(pBorn, ((e[0] * np.exp(-0.73 * e[2] * dist[i][j])/e[1])/(4 * np.pi * e[0] * dist[i][j]) - (e[2]/(2 * e[1] * (1 + (e[2]/2) * dist[i][j])))))
        pCor2 = np.append(pBorn,(1.0/(4* np.pi * e[0] * aa) * (((2.0 * (1.0 + u) * e[0] - (2 +2 * u + u**2)* e[1])/((1+u)* e[0] + (2 + 2*u + u**2)*e[1]))-(gamma+del1/(1+sig)) * (dist[i][j]/distInd[i][j]))))
        
potRF = potRF.reshape((len(Y), len(Y)))
potC = potC.reshape((len(Y), len(Y)))
pBorn = pBorn.reshape((len(Y), len(Y)))
pCor2 = pBorn.reshape((len(Y), len(Y)))

energy = []
for i in range(len(Y)):
    for j in range(len(Y)):
        energy = np.append(energy,(56 * (np.arcsin(np.power((1-cc**2/aa**2),0.5))-1) * (potC[i][j] + potRF[i][j] + pBorn[i][j] + pCor2[i][j])))
            
energy = energy.reshape((len(Y), len(Y)))
energy[np.isinf(energy)]= 0
energy[np.isnan(energy)]= 0 

Q = []
Pk=np.asarray([i[8] for i in S])
SA=np.asarray([i[9] for i in S])
Q=np.asarray([i[10] for i in S])
Restype=np.char.replace(np.char.replace(np.char.replace(np.char.replace(np.char.replace(Restype, 'HIS','H'), 'ASP','D'), 'ARG','R'), 'GLU','E'), 'LYS','K')
distSA = lambda u,v: (u+v)*0.5
SAij = 1.0-np.asarray([[distSA(u, v) for v in SA] for u in SA])

E_out = np.vstack([np.vstack([Q, -2.479*1.602*1.602*energy]), Pk])
np.savetxt('E.dat',E_out) 


if arguments.arg_s == 'MC': 
   
    print(u"\U0001F63A", "### TKSA - MC ###", u"\U0001F63A")
    start = time.time()        
    p = subprocess.Popen([r"c++","./src/tksamc.c",'-lm','-O3','-o','tksamc.exe'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
    p.communicate()
    p = subprocess.Popen(["./tksamc.exe",np.str(pH),np.str(T)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
    i=0
    j=1
    while p.poll() is None:
        sys.stdout.write('\r')
        sys.stdout.write("TKSA MC is running, please wait - [%20s%-20s]" % ('='*i,'='*i))
        sys.stdout.write(u"\U0001F63A")
        sys.stdout.flush()
        if i>19:
            j=j+1
        if j%2 == 0:
            i=i-1
        if j%2 == 1:
            i=i+1
        if i == 0:
            j=1
    sys.stdout.flush()
    time.sleep(0.1)
    
output,err = p.communicate()
print(output)
print(err)
end = time.time()
elapsed = end - start
print("Ran in %f sec" % elapsed)

 
if arguments.arg_plot == 'yes' and(arguments.arg_s =='MC'): 
    try: 
        file_plot = open("out.dat", 'r')
    except (IOError) as errno:
        print ('I/O error - ** Output file with issues - out.dat **. %s' % errno)
        sys.exit()

plot_data=[]
for line3 in file_plot: ## Plotting
    list3 = line3.split()
    plot_data.append(list3)

Restype=np.char.replace(np.char.replace(["%s%02d" % t for t in zip(Restype,np.asarray([i[1] for i in S]))],'C_T'+np.str(S[-1][1]),'CTR'),'N_T0'+np.str(S[0][1]),'NTR')
S=np.hstack((S,plot_data))

plot_data=list(map(float, np.asarray(plot_data).flatten()))
print("Total dG Energy: ",np.sum(np.asarray(plot_data)))
print(np.asarray(plot_data).reshape(len(Y),1))
width=1.0
x_pos = np.arange(len(total_charged_residues))+width/2.0
fig = plt.figure()
ax = fig.add_subplot(111)
colors = []
for position, value in enumerate(plot_data):
    if value > 0 and SA[position] > 0.5:
        colors.append('r')
    else:
        colors.append('b')
ax.bar(x_pos, plot_data,width=width,color=colors,linewidth=2,edgecolor='k')
ax.tick_params(axis='both',direction='in',bottom=True, top=True,right=True,left=True, length=5, width=2, which='major',labelsize=13)  
plt.setp(ax.spines.values(), linewidth=2)
if np.size(total_charged_residues)>35:
    plt.xticks(x_pos,Restype,rotation=90,fontsize=8)
elif np.size(total_charged_residues) >= 15 and np.size(total_charged_residues) <= 35:
    plt.xticks(x_pos,Restype,rotation=90,fontsize=12)
else: 
    plt.xticks(x_pos,Restype,rotation=90,fontsize=15)
plt.xlim([0,np.size(x_pos)])
plt.ylabel(r'$\Delta G_{qq}$(kJ/mol)',fontsize=20)
plt.show()
fig.savefig('Fig_MC_'+ os.path.splitext(arguments.f.name)[0]+'_pH_'+str(pH)+'_T_'+str(T)+'.jpg', dpi = 300)

header='1-Name	2-Residue-index	3-Position	4-Atom	5-Atom-type	6-X	7-Y	8-Z	9-PKA	10-SASA	11-Charge	12-dG_Energy 13-Total_dG= '+str(np.sum(np.asarray(plot_data)))+''
np.savetxt('Output_MC_'+os.path.splitext(arguments.f.name)[0]+'_pH_'+str(pH)+'_T_'+str(T)+'.dat',S,fmt='%s', delimiter="	",header=str(header))
cmd2 = 'mv result.txt *.exe E.dat out.dat SASA* x.txt y.txt z.txt energiaxyza.txt energiaxyza_fixed.txt '+os.path.splitext(arguments.f.name)[0]+'*.txt ./aux'
os.system(cmd2)
# print(u"\U0001F63A", "### Finished ###", u"\U0001F63A")




#def plot_ellipsoid_3d(ell, ax):
#    """Plot the 3-d Ellipsoid ell on the Axes3D ax."""

    # points on unit sphere
#    u = np.linspace(0.0, 2.0 * np.pi, 100)
#    v = np.linspace(0.0, np.pi, 100)
#    z = np.outer(np.cos(u), np.sin(v))
#    y = np.outer(np.sin(u), np.sin(v))
#    x = np.outer(np.ones_like(u), np.cos(v))

    # transform points to ellipsoid
#    for i in range(len(x)):
#        for j in range(len(x)):
#            x[i,j], y[i,j], z[i,j] = ell.ctr + np.dot(ell.axes,
 #                                                     [x[i,j],y[i,j],z[i,j]])

 #   ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color='#2980b9', alpha=0.4)


# Generate points within a torus




# Generate samples drawn from a single ellipsoid
#npoints = len(Y)
#points = XYZ

#pointvol = ell_gen.vol / npoints

# Find bounding ellipsoid(s)
#ells = nestle.bounding_ellipsoids(points, pointvol)


# plot
#fig = plt.figure(figsize=(10., 10.))
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(points[:, 0], points[:, 1], points[:, 2], c='k', marker='o')
#for ell in ells:
#    plot_ellipsoid_3d(ell, ax)

#ax.set_xlabel('X Label')
#ax.set_ylabel('Y Label')
#ax.set_zlabel('Z Label')

#ax.set_xlim((-np.max(dist_origin),np.max(dist_origin)))
#ax.set_ylim((-np.max(dist_origin),np.max(dist_origin)))
#ax.set_zlim((-np.max(dist_origin),np.max(dist_origin)))
#fig.tight_layout()
#plt.show()



