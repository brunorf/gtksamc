#!/usr/bin/env python
#coding: utf8 

__description__ = \
"""
Charge-charge energy calculation in python
"""

__author__ = "Dr. Vinícius Contessoto"
__date__ = "21/12/2016"

################################################################
# TKSA em python - faz tudo sem precisar do bash
#
# Try to rewrite in python way by Vinícius Contessoto
# start 30/09/2015
#
# Nesta versão, faz a leitura em python mas roda em c++
#
# python pulo_do_gato.py PDB.pdb
# 
# Need to install some bib from python
################################################################

import sys
import numpy as np
import scipy as sc
import subprocess
import os
import argparse
import time
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Plot dGqq vs pH in python')
parser.add_argument('-f', metavar='input-file-Gqq-ph',help='insert dGqq-ph file',type=argparse.FileType('rt'))

try:
    arguments = parser.parse_args()
    #print '################################################'
    #print u"\U0001F63A", "### Plot Pulo do Gato iniciado ###", u"\U0001F63A"
    #print 'Input file:', arguments.f.name
    
except IOError, msg:
    parser.error(str(msg))                    
                  

def main():

   file_plot = arguments.f     
   
   #try: 
   #  file_plot = open(arguments.f, 'r')
   #except (IOError) as errno:
   #  print ('I/O error - ** Arquivo out.dat com problemas **. %s' % errno)
   #  sys.exit()
   
   plot_data=[]
   for line in file_plot: ## Calculando a SASA
     list1 = line.split()
     plot_data.append(list1)
   
   #print plot_data
   #Restype=np.char.replace(np.char.replace(["%s%02d" % t for t in zip(Restype,np.asarray([i[1] for i in S]))],'C_T'+np.str(S[-1][1]),'CTR'),'N_T0'+np.str(S[0][1]),'NTR')
   #S=np.hstack((S,plot_data))
   #
   
   #plot_data=list(map(float, np.asarray(plot_data).flatten()))
   #print plot_data
   #print "Total dG Energy: ",np.sum(np.asarray(plot_data))
   #x_pos = np.arange(len(total_charged_residues))
   fig = plt.figure()
   ax = fig.add_subplot(111)
   #width=1.0
   #colors = []
   #for position, value in enumerate(plot_data):
   #  if value > 0 and SA[position] > 0.5:
   #     colors.append('r')
   #  else:
   #     colors.append('b')
   #ax.scatter(plot_data,linewidth=2)
   #ax.tick_params('both', length=5, width=2, which='major',labelsize=13)   
   #plt.setp(ax.spines.values(), linewidth=2)
   #plt.xticks(x_pos+width/2.0,Restype,rotation=90,fontsize=15)
   #plt.xlim([0,np.size(x_pos)])
   ##plt.ylim([-12.0,4.0])
   #plt.ylabel(r'$\Delta G_{qq}$(kJ/mol)',fontsize=20)
   ##plt.title('Protein ' + os.path.splitext(arguments.f.name)[0])

   X=np.asarray([i[0] for i in plot_data]).astype(np.float)
   Y=np.asarray([i[1] for i in plot_data]).astype(np.float)

   plt.plot(X,Y,'-o',markersize=10,linewidth=2.5)
   ax.tick_params('both', length=5, width=2.5, which='major',labelsize=15)
   plt.setp(ax.spines.values(), linewidth=2)
   plt.xticks(np.arange(min(X)-1.0, max(X)+2.0, 1.0))

   plt.xlabel(r'pH',fontsize=20)
   plt.ylabel(r'$\Delta G_{elec}$(kJ/mol)',fontsize=20)
   #plt.ylim([-65,-51])
   #plt.show()
   fig.savefig('Fig_Gqq-pH_'+ os.path.splitext(arguments.f.name)[0]+'.jpg', dpi = 300)
   #
   #header='1-Name	2-Residue-index	3-Position	4-Atom	5-Atom-type	6-X	7-Y	8-Z	9-PKA	10-SASA	11-Charge	12-dG_Energy 13-Total_dG= '+str(np.sum(np.asarray(plot_data)))+''
   #Y=plot_data[:,2]
   #Y=plot_data[:,2]
   #np.savetxt('Output_EX_'+os.path.splitext(arguments.f.name)[0]+'_pH_'+str(pH)+'_T_'+str(T)+'.dat',S,fmt='%s', delimiter="	",header=str(header))
   
   #print u"\U0001F63A", "### Finished ###", u"\U0001F63A"
      
if __name__ == "__main__": main()
