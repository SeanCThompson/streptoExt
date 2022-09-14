#!/usr/bin/python3

'''
Makes plots of the metabolic switch transition probabilities.
Puts on the x axis a parameter of your choice. 
The parameter you want to use must be in the file name, and you have to specify its position or the subtring

Usage: 
./plot_switching_params.py [file]

'''

import matplotlib as mpl
mpl.use('Qt5Agg')

import sys,math,os,string,random
from subprocess import Popen, PIPE
from collections import Counter

#from PIL import Image
# import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import numpy as np

ltime=[]
l_p_g_a=[]
l_p_a_g=[]
l_avg_p_g_a=[]
l_avg_p_a_g=[]
l_std_p_g_a=[]
l_std2_p_g_a=[]
l_std_p_a_g=[]
l_std2_p_a_g=[]
l_p_g_a_stress=[]
l_p_a_g_stress=[]
l_avg_p_g_a_stress=[]
l_avg_p_a_g_stress=[]
l_std_p_g_a_stress=[]
l_std2_p_g_a_stress=[]
l_std_p_a_g_stress=[]
l_std2_p_a_g_stress=[]
t_init=0

genome_pos_in_file = 4

filename = sys.argv[1]
maxtime =float("inf") 

try :
    maxtime = int(sys.argv[2])
except:
    pass
time = int( (Popen(["head","-n1",filename], stdout=PIPE).stdout.read()).split()[0])
print( "Initial time =",time)
with open(filename,"r") as fin:
    for line in fin:
        line_space =line.split()
        if len(line_space)< 4: continue
        if "n" in line_space: continue

        line_comma = line.split(", ")
        

        try:
            trans_prob=line_comma[2]
            trans_prob = trans_prob.split()
            l_p_g_a.append(float(trans_prob[0].replace(",","")))
            l_p_a_g.append(float(trans_prob[1].replace(",","")))
                        
            trans_prob=line_comma[3]
            trans_prob = trans_prob.split()
            l_p_g_a_stress.append(float(trans_prob[0].replace(",","")))
            l_p_a_g_stress.append(float(trans_prob[1].replace(",","")))
        except:
            pass

        timenow=int(line_space[0])
        if timenow != time:
            if time>maxtime: break
            


            ltime.append(time)
            if time%50000==0: print(time)
            time = timenow
            l_avg_p_g_a.append(np.median(l_p_g_a))
            l_std_p_g_a.append(np.quantile(l_p_g_a,0.25))
            l_std2_p_g_a.append(np.quantile(l_p_g_a,0.75))
            l_avg_p_a_g.append(np.median(l_p_a_g))
            l_std_p_a_g.append(np.quantile(l_p_a_g,0.25))
            l_std2_p_a_g.append(np.quantile(l_p_a_g,0.75))
            l_p_g_a=[]
            l_p_a_g=[]

            l_avg_p_g_a_stress.append(np.median(l_p_g_a_stress))
            l_std_p_g_a_stress.append(np.quantile(l_p_g_a_stress,0.25))
            l_std2_p_g_a_stress.append(np.quantile(l_p_g_a_stress,0.75))
            l_avg_p_a_g_stress.append(np.median(l_p_a_g_stress))
            l_std_p_a_g_stress.append(np.quantile(l_p_a_g_stress,0.25))
            l_std2_p_a_g_stress.append(np.quantile(l_p_a_g_stress,0.75))
            l_p_g_a_stress=[]
            l_p_a_g_stress=[]

fig = plt.figure(filename)
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2,sharex=ax1, sharey=ax1)
#fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)
#fig(filename)

ax1.plot(ltime,l_avg_p_g_a,label = "median P(g to a) +/- 25\%")
ax1.plot(ltime,l_avg_p_a_g,label = "median P(a to g) +/- 25\%")

ax1.fill_between(ltime, l_std_p_g_a,l_std2_p_g_a , alpha=0.5)
ax1.fill_between(ltime, l_std_p_a_g,l_std2_p_a_g,alpha=0.5)

ax1.set_xlim(xmin=0)
ax1.set_xlabel("Time steps")
ax1.set_ylabel("Trans Prob")

logscale=True
if logscale:
    #plt.yscale('log')
    ax1.set_ylim(ymin=0)
else:
    ax1.set_ylim(ymin=0)
    
ax1.legend()
ax1.set_title("Nominal")

ax2.plot(ltime,l_avg_p_g_a_stress,label = "median P(g to a) +/- 25\%")
ax2.plot(ltime,l_avg_p_a_g_stress,label = "median P(a to g) +/- 25\%")

ax2.fill_between(ltime, l_std_p_g_a_stress,l_std2_p_g_a_stress, alpha=0.5)
ax2.fill_between(ltime, l_std_p_a_g_stress,l_std2_p_a_g_stress,alpha=0.5)

ax2.set_xlim(xmin=0)
ax2.set_xlabel("Time steps")
ax2.set_ylabel("Trans Prob")

logscale=True
if logscale:
    #plt.yscale('log')
    ax2.set_ylim(ymin=0)
else:
    ax2.set_ylim(ymin=0)
    
ax2.legend()
ax2.set_title("Stressed")

plt.show()


