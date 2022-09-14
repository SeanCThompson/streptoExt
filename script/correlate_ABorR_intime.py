#!/usr/bin/python3

'''
Makes plots of genome and genes distribution.
Positions must be specified in some way.. for now hard coded
-also added antib production
'''
import matplotlib as mpl
mpl.use('Qt5Agg')

import sys,math,os,string,random,copy
from subprocess import Popen, PIPE
from collections import Counter

#from PIL import Image
# mpl.use('GTK3Agg')
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import numpy as np
import seaborn as sns


def Gene2Number(gene):
    if gene=='F':return 0
    if gene=='A':return 1
    if gene=='B':return 2

# dg={}
# dg_nominal={}
# dg_stress={}
count=0
broken_repl=0
filename = sys.argv[1]
# tot_people_alive_pergnm={}

mintime=2000
maxtime = 2500 # 1000#2500

def extractWindow(mintime, maxtime, filename):

    dg={}
    dg_nominal={}
    dg_stress={}
    line_number=0
    stressCount=0
    with open(filename,"r") as fin:
        
        for line in fin:
            line_number+=1
            line = line.split()
            time = int(line[0])
            if int(line[-1])==1: continue
            if time > maxtime : break
            
            # E = extinction
            if "E:" in line[1]:
                #we be done early!
                for key in dg:
                    if dg[key]["death"]==maxtime:
                        dg[key]["death"]=time+1 # if someone is already died we don't artificially extend its life time
                break
            tag=line[2]
            if tag not in dg:
                #set line 2
                # I = initial 
                if "I:" in line[1]:
                    try:
                        dg[tag]={"gnm":line[3], "birth":time, "death":maxtime, "A":0, "R":0, "offs_tag":[], "parent":"n"} #set death to 2500, this will be true unless death is recorded in log
                        dg_nominal[tag]={"gnm":line[3], "birth":time, "death":maxtime, "A":0, "R":0, "offs_tag":[], "parent":"n"} # The initial state doesn't have a stress state that is meaningful
                        dg_stress[tag]={"gnm":line[3], "birth":time, "death":maxtime, "A":0, "R":0, "offs_tag":[], "parent":"n"} # The initial state doesn't have a stress state that is meaningful
                    except:
                        print("Error, here is line: ", line)
                        print("line_number = ", line_number)
                        # sys.exit(1)
                        pass
                else:
                    if len(line)>3:
                        print("Error, this is not an initial guy.")
                        print(line)
                        #sys.exit(1)
                    else:
                        continue # it's a guy with zero length genome
            else:
                #update life
                # R = replication
                if "R:" in line[1]:
                    #it replicated
                    if time >= mintime: 
                        dg[tag]["R"]+=1
                    if len(line)>5:
                        offspr_tag = line[4]
                        dg[tag]["offs_tag"].append(offspr_tag)
                        if line[6] == "1":
                            if time==maxtime:
                                stressCount+=1
                            dg_stress[tag]["offs_tag"].append(offspr_tag)
                            if time >= mintime:
                                dg_stress[tag]["R"]+=1
                                dg_stress[tag]["gnm"]=line[3]
                        else:
                            dg_nominal[tag]["offs_tag"].append(offspr_tag)
                            if time >= mintime:
                                dg_nominal[tag]["R"]+=1
                        if offspr_tag in dg:
                            print("Error, this offspring tag is already in dg")
                        else:
                            if len(line)>5:
                                dg[offspr_tag]={"gnm":line[5], "birth":time, "death":maxtime, "A":0, "R":0, "offs_tag":[], "parent":tag} #set death to maxtime, this will be true unless death is recorded in log
                                dg_stress[offspr_tag]={"gnm":line[5], "birth":time, "death":maxtime, "A":0, "R":0, "offs_tag":[], "parent":tag} #set death to maxtime, this will be true unless death is recorded in log
                                dg_nominal[offspr_tag]={"gnm":line[5], "birth":time, "death":maxtime, "A":0, "R":0, "offs_tag":[], "parent":tag} #set death to maxtime, this will be true unless death is recorded in log
                # A = antibiotic production
                elif "A:" in line[1]:
                    #it made an antibiotic
                    if time >= mintime:
                        dg[tag]["A"]+=1
                    if line[5]=="1":
                        if time >= mintime: 
                            if time==maxtime:
                                stressCount+=1
                            dg_stress[tag]["A"]+=1
                            dg_stress[tag]["gnm"]=line[3]
                    else:
                        if time >= mintime: 
                            dg_nominal[tag]["A"]+=1
                # D = death
                elif "D:" in line[1]:
                    if tag in dg:
                        dg[tag]["death"]=time
                        if line[4]=="1":
                            if time==maxtime:
                                stressCount+=1
                            dg_stress[tag]["death"]=time
                            dg_stress[tag]["gnm"]=line[3]
                        else:
                            dg_nominal[tag]["death"]=time

    return dg, dg_nominal, dg_stress, stressCount
        
dg, dg_nominal, dg_stress, stressCount = extractWindow(1000, 1100, filename)

print("Done, reading data. In total so many people lived:", len(dg))




def set_size(w,h, ax=None):
    """ w, h: width, height in inches """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)



print ("Hello 1")
# max_nr_F=0
# for key in dg:
#     max_nr_F = max(dg[key]["gnm"].count("F") , max_nr_F)

max_nr_Genes=0
for key in dg:
    max_nr_Genes = max(len(dg[key]["gnm"]) , max_nr_Genes)

def extractLabor(dg, max_nr_Genes):
    lab=[0 for _ in range(max_nr_Genes+1)]
    lrepl=[0 for _ in range(max_nr_Genes+1)]
    for key in dg:
        max_nr_Genes = len(dg[key]["gnm"])
        lab[ max_nr_Genes ]+= dg[key]["A"] 
        lrepl[ max_nr_Genes ]+= dg[key]["R"]
        
    return lab, lrepl


from operator import add
ab_mat=np.zeros(shape=(25,max_nr_Genes+1), dtype=float)
repl_mat=np.zeros(shape=(25,max_nr_Genes+1), dtype=float)
l_stress_ratio=[]
row=0
for time in range(0, 2500, 100):
    dg, dg_nominal, dg_stress, stressCount = extractWindow(time, time+100, filename)
    #print(stressCount/len(dg_nominal))
    l_stress_ratio.append(stressCount/len(dg_nominal))
    lab_nominal, lrepl_nominal =extractLabor(dg_nominal, max_nr_Genes)
    lab_stress, lrepl_stress =extractLabor(dg_stress, max_nr_Genes)



    lab = list( map(add, lab_nominal, lab_stress) )
    lrepl = list( map(add, lrepl_nominal, lrepl_stress) )

    ab_mat[row]=[x/float(sum(lab)+1) for x in lab]
    repl_mat[row]=[x/float(sum(lrepl)+1) for x in lrepl]
    ab_mat[row]=[x for x in lab]
    repl_mat[row]=[x for x in lrepl]
    #ab_mat[row]=lab
    row+=1

my_cmapF = sns.color_palette("gist_earth_r", as_cmap=True)
# F genes
#X, Y = np.meshgrid( range(maxlen+1) , np.linspace(0,bin_gnm,num=bin_gnm))


lab, lrepl =extractLabor(dg, max_nr_Genes)
#print(l_stress_ratio)
#print(ab_mat)
print(ab_mat[-1])
#print(lab_stress)
#print(lab_nominal)

# hist_df, bin_edgesdf = np.histogram( ldeltaF, bins = np.linspace(0.,1.05,num=21) )
# hist_df = [x/float(len(dg)) for x in hist_df]


print("Plotting")

# fig, ax = plt.subplots(2,2)
fig = plt.figure()
fig.canvas.set_window_title(sys.argv[1])
gs = gridspec.GridSpec(6, 2)
ax=[]
ax.append([fig.add_subplot(gs[0:2, 0:2]),fig.add_subplot(gs[2:4,0:2]), fig.add_subplot(gs[4:6,0:2])])

#print(ax)

sum_lab=sum(lab)
sum_repl=sum(lrepl)

sum_lab_nominal=sum(lab_nominal)
sum_repl_nominal=sum(lrepl_nominal)

sum_lab_stress=sum(lab_stress)
sum_repl_stress=sum(lrepl_stress)



# ax[0][0].plot(range(max_nr_F+1),[x/float(sum_lab) for x in lab_nominal],label="AB produced (nominal)", color='blue')
# ax[0][0].plot(range(max_nr_F+1),[x/float(sum_repl) for x in lrepl_nominal],label="Replication (nominal)", color='orange')

# ax[0][0].plot(range(max_nr_F+1),[x/float(sum_lab) for x in lab_stress],label="AB produced (stress)", linestyle = 'dotted', color='blue')
# ax[0][0].plot(range(max_nr_F+1),[x/float(sum_repl) for x in lrepl_stress],label="Replication (stress)", linestyle = 'dotted', color='orange')




X, Y = np.meshgrid( range(max_nr_Genes+1) , range(0,25))
ab_2dplot = ax[0][0].pcolormesh(Y,X, repl_mat, cmap = my_cmapF )
cb1=fig.colorbar(ab_2dplot, ax=ax[0][0], orientation='horizontal', anchor=(0.0, 0.25), shrink=0.15, fraction=0.075)
cb1.ax.tick_params(labelsize='small', length=2)
ax[0][0].set_xlim(xmin=0, xmax=24)
ax[0][0].set_ylabel("Genome length")
ax[0][0].set_title("Reproduction events")
spine=ax[0][0].spines["top"]
spine.set_visible(False)
spine=ax[0][0].spines["right"]
spine.set_visible(False)

F_2dplot = ax[0][1].pcolormesh(Y,X, ab_mat, cmap = my_cmapF )
cb2=fig.colorbar(F_2dplot, ax=ax[0][1], orientation='horizontal', anchor=(0.0, 0.25), shrink=0.15, fraction=0.075)
cb2.ax.tick_params(labelsize='small', length=2)
ax[0][1].set_xlim(xmin=0, xmax=24)
ax[0][1].set_ylabel("Genome length")
ax[0][1].set_title("Antibiotic production events")
spine=ax[0][1].spines["top"]
spine.set_visible(False)
spine=ax[0][1].spines["right"]
spine.set_visible(False)

ax[0][2].plot(range(0, 25),[x*100 for x in l_stress_ratio])
ax[0][2].set_xlabel("Time (100 steps)")
ax[0][2].set_ylabel("% of colony in stressed state")
ax[0][2].set_xlim(xmin=0, xmax=24)
ax[0][2].set_title("Stressed percentage of colony")

plt.show()

print(sum(ab_mat[-1][0:(len(ab_mat[-1])-1)]))
print((sum(ab_mat[-1][0:(len(ab_mat[-1])-1)])/sum(ab_mat[-1]))*100)

plt.plot(ab_mat[-1][0:(len(ab_mat[-1])-1)])
plt.ylabel('Antibiotic producton events during season')
plt.xlabel('Genome length')
plt.show()
