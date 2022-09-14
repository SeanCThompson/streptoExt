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


def Gene2Number(gene):
    if gene=='F':return 0
    if gene=='A':return 1
    if gene=='B':return 2

dg={}
dg_nominal={}
dg_stress={}
count=0
broken_repl=0
filename = sys.argv[1]
# tot_people_alive_pergnm={}
line_number=0
mintime=2000
maxtime = 2500 # 1000#2500

with open(filename,"r") as fin:
    
    for line in fin:
        line_number+=1
        line = line.split()
        time = int(line[0])
        if int(line[-1])==1: continue
        if time >= maxtime : break
        
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
                dg[tag]["R"]+=1
                if len(line)>5:
                    offspr_tag = line[4]
                    dg[tag]["offs_tag"].append(offspr_tag)
                    if line[6] == "1":
                        
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
                dg[tag]["A"]+=1
                if line[5]=="1":
                    if time >= mintime: 
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
                        dg_stress[tag]["death"]=time
                        dg_stress[tag]["gnm"]=line[3]
                    else:
                        dg_nominal[tag]["death"]=time
        

print("Done, reading data. In total so many people lived:", len(dg))

print("Successful genomes - in order:")
l_nr_offs=[]
d_genomes={} #dictionary with genomes a keys; values: total abundance, earliest found
d_gnm_and_offs={} #dict with parent offspring relation. key: genome; values: descending genomes
for key in dg:
    gnm = dg[key]["gnm"]
    if gnm in d_genomes:
        d_genomes[gnm]["howmany_offs"]+=len(dg[key]["offs_tag"])
        d_genomes[gnm]["earliest"] = min( dg[key]["birth"] , d_genomes[gnm]["earliest"] )
        d_genomes[gnm]["howmany"]+=1
        d_genomes[gnm]["ab_prod"]+=dg[key]["A"]
    else: d_genomes[gnm]={"ab_prod": dg[key]["A"],  "howmany": 1 , "howmany_offs":len(dg[key]["offs_tag"]) , "earliest": dg[key]["birth"] }

    if gnm in d_gnm_and_offs:
        d_gnm_and_offs[gnm].extend( [dg[ot]["gnm"] for ot in dg[key]["offs_tag"]] )
    else:
        d_gnm_and_offs[gnm] = [ dg[ot]["gnm"] for ot in dg[key]["offs_tag"] ]




#just mutational relations
for gnm in d_gnm_and_offs:
    d_gnm_and_offs[gnm]=set(d_gnm_and_offs)



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
max_nr_F=0
for key in dg:
    max_nr_F = max(dg[key]["gnm"].count("F") , max_nr_F)

print ("Hello 3")
#first thing to check: Who makes AB and who reproduces:
lab=[0 for _ in range(max_nr_F+1)]
lrepl=[0 for _ in range(max_nr_F+1)]
# print(lab)
ldeltaF=[]
ldeltaF2=[]
# d_genomes={}
for key in dg:
    me_nr_F = dg[key]["gnm"].count("F")
    lab[ me_nr_F ]+= dg[key]["A"] 
    lrepl[ me_nr_F ]+= dg[key]["R"]
    
    # me_genome = dg[key]["gnm"]
    # if me_genome in d_genomes:
    #     d_genomes[me_genome]["howmany"]+=1
    # else:
    #     d_genomes[me_genome]={"howmany":1, "offs_gnm":[]}
    
    if len(dg[key]["offs_tag"])>0:
        this_deltaf=[]
        for offs_tag in dg[key]["offs_tag"]:
            # print( dg[key]["gnm"] )
            # print( dg[offs_tag]["gnm"] )
            deltaFfract = ( dg[key]["gnm"].count("F") - dg[offs_tag]["gnm"].count("F") )/(float( dg[key]["gnm"].count("F") ))
            ldeltaF.append( deltaFfract )
            # this_deltaf.append( deltaFfract )
            
            # d_genomes[me_genome]["offs_gnm"].append( dg[offs_tag]["gnm"] )
            # if dg[offs_tag]["gnm"] in d_genomes:
            #     d_genomes[ dg[offs_tag]["gnm"] ]["howmany"]+=1
            # else:
            #     d_genomes[ dg[offs_tag]["gnm"] ]={"howmany":1, "offs_gnm":[]}
        # len_thisdeltaf = float(len(this_deltaf))
        # hist_thisdeltaf,bins = np.histogram( this_deltaf, bins = np.linspace(0.,1.1,num=11) )
        # ldeltaF2.append([x/len_thisdeltaf for x in hist_thisdeltaf])

def extractLabor(dg, max_nr_F):
    lab=[0 for _ in range(max_nr_F+1)]
    lrepl=[0 for _ in range(max_nr_F+1)]
    ldeltaF=[]
    for key in dg:
        me_nr_F = dg[key]["gnm"].count("F")
        lab[ me_nr_F ]+= dg[key]["A"] 
        lrepl[ me_nr_F ]+= dg[key]["R"]
        
        
        # if len(dg[key]["offs_tag"])>0:
        #     this_deltaf=[]
        #     for offs_tag in dg[key]["offs_tag"]:

        #         #deltaFfract = ( dg[key]["gnm"].count("F") - dg[offs_tag]["gnm"].count("F") )/(float( dg[key]["gnm"].count("F") ))
        #         #ldeltaF.append( deltaFfract )
    return lab, lrepl

lab_nominal, lrepl_nominal =extractLabor(dg_nominal, max_nr_F)
lab_stress, lrepl_stress =extractLabor(dg_stress, max_nr_F)

print(lrepl_stress)


hist_df, bin_edgesdf = np.histogram( ldeltaF, bins = np.linspace(0.,1.05,num=21) )
hist_df = [x/float(len(dg)) for x in hist_df]


print("Plotting")

# fig, ax = plt.subplots(2,2)
fig = plt.figure()
fig.canvas.set_window_title(sys.argv[1])
gs = gridspec.GridSpec(4, 4)
ax=[]
ax.append([fig.add_subplot(gs[0:2, 0:2]), fig.add_subplot(gs[0:2,2:4])] )
ax.append([fig.add_subplot(gs[2:4, 0:2]),fig.add_subplot(gs[2:4,2]),fig.add_subplot(gs[2:4,3]) ])

print(ax)

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


ax[0][0].plot(range(max_nr_F+1),[x for x in lab_nominal],label="AB produced (nominal)", color='blue')
#ax[0][0].plot(range(max_nr_F+1),[x/float(sum_repl) for x in lrepl_nominal],label="Replication (nominal)", color='orange')

ax[0][0].plot(range(max_nr_F+1),[x for x in lab_stress],label="AB produced (stress)", linestyle = 'dotted', color='blue')
#ax[0][0].plot(range(max_nr_F+1),[x/float(sum_repl) for x in lrepl_stress],label="Replication (stress)", linestyle = 'dotted', color='orange')

#ax[1][0].plot(range(max_nr_F+1),[x/float(sum_lab) for x in lab_nominal],label="AB produced (nominal)", color='blue')
ax[1][0].plot(range(max_nr_F+1),[x for x in lrepl_nominal],label="Replication (nominal)", color='orange')

#ax[1][0].plot(range(max_nr_F+1),[x/float(sum_lab) for x in lab_stress],label="AB produced (stress)", linestyle = 'dotted', color='blue')
ax[1][0].plot(range(max_nr_F+1),[x for x in lrepl_stress],label="Replication (stress)", linestyle = 'dotted', color='orange')




ax[0][0].set_xlabel("Nr. of Growth-promoting genes")
ax[0][0].legend()

ax[0][1].plot(bin_edgesdf[:-1],hist_df, label="density")
ax[0][1].set_xlabel("Fraction of gr-promoting genes after mutations")
# ax[1].set_yscale('log')
# print("Here1")
# ax[1].plot(bin_edgesdf[:-1],lmedian_deltaf2, label="median density 2")
# ax[1].plot(bin_edgesdf[:-1],lmean_deltaf2, label="mean density 2")
# ax[1].plot(bin_edgesdf[:-1],l25per_deltaf2, label="25perc density 2")
# ax[1].plot(bin_edgesdf[:-1],l75per_deltaf2, label="75 per 2")
# print("Here2")
ax[0][1].legend()

ax[0][1].fill_between(bin_edgesdf[:-1],hist_df,[0. for x in bin_edgesdf[:-1]], alpha=0.3)
ax[0][1].plot( [ 0 , 1 ] , [0,0], lw=0.5, c='red', linestyle = 'dotted' ) #a zero line to facilitate comparison
ax[0][1].set_xlim([0.,1.])


# print("Making Hist")
# H, xedges, yedges = np.histogram2d(ldeltaf_po, ldeltaab, bins=(range(20), range(-10,10,1)) )
# print("Transposing")
# H = H.T
# print("Meshgrid")
# xx, yy = np.meshgrid(xedges, yedges)
# # mask some 'bad' data, in your case you would have: data == 0
# H = np.ma.masked_where(H < 0.0000001, H)
# # cmap = plt.cm.magma_r
# cmap = copy.copy(mpl.cm.get_cmap("magma_r"))
# cmap.set_bad(color='white')
# print("plot")
# pcm = ax[2].pcolormesh(xx, yy, H,cmap=cmap)

# ax[2].plot(range(1+max(ldeltaf_po)), lavr_deltaab )
# ax[2].fill_between(range(1+max(ldeltaf_po)), [x+y for x,y in zip(lavr_deltaab,lstd_deltaab)],[x-y for x,y in zip(lavr_deltaab,lstd_deltaab)],alpha=0.3 )
# ax[2].plot(range(1+max(ldeltaf_po)), [x-y for x,y in zip(lavr_deltaab,lstd_deltaab)] )
# ax[1][0].plot(range(1+max(ldeltaf_po)), lmedian_deltaab, label=r'$\Delta$'+" antibiotic production")
# ax[1][0].fill_between(range(1+max(ldeltaf_po)), l25perc_deltaab, l75perc_deltaab, alpha=0.3 )
# ax[1][0].plot( [ 0 , max(ldeltaf_po) ] , [0,0], lw=0.5, c='red', linestyle = 'dotted' ) #a zero line to facilitate comparison
# ax[1][0].set_xlabel(r'$\Delta$'+" nr. growth-promoting genes")
# ax[1][0].legend()
# fig.colorbar(pcm)
# c = ax.imshow(H,interpolation ='none', origin ='lower') 
# plt.tight_layout()
#print("l_repr_ab_for_lollyplot",l_repr_ab_for_lollyplot)
#makes sense to make lollyplots only if we got founders, otherwise it's a bit messy
if len(sys.argv)>2:
    #lollyplot
    # l_repr_ab_for_lollyplot= sorted( l_repr_ab_for_lollyplot, key=lambda x: (x[1],-x[0]),reverse=True ) 
    l_repr_ab_for_lollyplot=l_repr_ab_for_lollyplot[::-1]
    my_range = range( len(l_repr_ab_for_lollyplot) )
    I_want_log10 = False
    I_want_fractions = False
    I_want_absolute_numbers=False
    I_am_hacking_stuff = False
    I_want_scatterplot = False
    I_want_separate_plots = True
    if I_want_log10:
        ax[1][1].hlines(y=my_range, xmin=[ -math.log10(x[0]) if x[0]>0 else 0. for x in l_repr_ab_for_lollyplot ], 
                                                    xmax=[  math.log10(x[1]) if x[1]>0 else 0. for x in l_repr_ab_for_lollyplot ], 
                                                    color='grey', alpha=.8)
        ax[1][1].plot([0,0],[min(my_range),max(my_range)], lw=0.5, c='red', linestyle = 'dotted')
        ax[1][1].scatter([ -math.log10(x[0]) if x[0]>0 else 0. for x in l_repr_ab_for_lollyplot ], my_range, color='blue', alpha=1., label='repl')
        ax[1][1].scatter([  math.log10(x[1]) if x[1]>0 else 0. for x in l_repr_ab_for_lollyplot ] , my_range, color='orange', alpha=1. , label='AB')
        ax[1][1].legend()
    elif I_want_absolute_numbers:
        ax[1][1].hlines(y=my_range, xmin=[ -x[0] for x in l_repr_ab_for_lollyplot ], 
                                    xmax=[  x[1] for x in l_repr_ab_for_lollyplot ], 
                                    color='grey', alpha=.8)
        ax[1][1].plot([0,0],[min(my_range),max(my_range)], lw=0.5, c='red', linestyle = 'dotted')
        ax[1][1].scatter([ -x[0] for x in l_repr_ab_for_lollyplot ], my_range, color='blue', alpha=1., label='repl')
        ax[1][1].scatter([  x[1] for x in l_repr_ab_for_lollyplot ], my_range, color='orange', alpha=1. , label='AB')
        ax[1][1].legend()
    elif I_want_fractions:
        sum_l_repr_ab_for_lollyplot_0 = sum(x[0] for x in l_repr_ab_for_lollyplot)
        sum_l_repr_ab_for_lollyplot_1 = sum(x[1] for x in l_repr_ab_for_lollyplot)
        l_repr_ab_for_lollyplot = [ [x[0]/float(sum_l_repr_ab_for_lollyplot_0), x[1]/float(sum_l_repr_ab_for_lollyplot_1) ] for x in l_repr_ab_for_lollyplot ]
        ax[1][1].hlines(y=my_range, xmin=[ -x[0] for x in l_repr_ab_for_lollyplot ], 
                                    xmax=[  x[1] for x in l_repr_ab_for_lollyplot ], 
                                    color='grey', alpha=.8)
        ax[1][1].plot([0,0],[min(my_range),max(my_range)], lw=0.5, c='red', linestyle = 'dotted')
        ax[1][1].scatter([ -x[0] for x in l_repr_ab_for_lollyplot ], my_range, color='blue', alpha=1., label='repl')
        ax[1][1].scatter([  x[1] for x in l_repr_ab_for_lollyplot ], my_range, color='orange', alpha=1. , label='AB')
        ax[1][1].legend()
    elif I_want_scatterplot:
        sum_l_repr_ab_for_lollyplot_0 = sum(x[0] for x in l_repr_ab_for_lollyplot)
        sum_l_repr_ab_for_lollyplot_1 = sum(x[1] for x in l_repr_ab_for_lollyplot)
        l_repr_ab_for_lollyplot = [ [x[0]/float(sum_l_repr_ab_for_lollyplot_0), x[1]/float(sum_l_repr_ab_for_lollyplot_1),x[2] ] for x in l_repr_ab_for_lollyplot ]
        
        # ax[1][1].scatter([ x[1] if x[1]>0 else 0. for x in l_repr_ab_for_lollyplot ], 
        #                 [ x[0] if x[0]>0 else 0. for x in l_repr_ab_for_lollyplot ] )
        # sorted_to_plot = l_repr_ab_for_lollyplot.sort(key=lambda x: x[1]) # sort on AB
        bla = ax[1][1].scatter([ x[0] for x in l_repr_ab_for_lollyplot ], 
                        [ x[1] for x in l_repr_ab_for_lollyplot ])
        print("l_repr_ab_for_lollyplot",l_repr_ab_for_lollyplot)
        for i,x in enumerate(l_repr_ab_for_lollyplot):
            print("x",x)
            ax[1][1].annotate(x[2], (x[0],x[1]),fontsize=2 ,rotation=45)
        # ax[1][1].invert_xaxis()
        # ax[1][1].invert_yaxis()
        ax[1][1].set_ylabel("AB")
        ax[1][1].set_xlabel("repl")
    if I_want_separate_plots:
        #we makes two separate plots in ax[1][1] and ax[1][2]
        sum_l_repr_ab_for_lollyplot_0 = sum(x[0] for x in l_repr_ab_for_lollyplot)
        sum_l_repr_ab_for_lollyplot_1 = sum(x[1] for x in l_repr_ab_for_lollyplot)
        l_repr_ab_for_lollyplot = [ [x[0]/float(sum_l_repr_ab_for_lollyplot_0), x[1]/float(sum_l_repr_ab_for_lollyplot_1) ] for x in l_repr_ab_for_lollyplot ]
        ax[1][1].hlines(y=my_range, xmin=[ 0 for x in l_repr_ab_for_lollyplot ], 
                                    xmax=[  x[0] for x in l_repr_ab_for_lollyplot ], 
                                    color='grey', alpha=.8)
        ax[1][2].hlines(y=my_range, xmin=[ 0 for x in l_repr_ab_for_lollyplot ], 
                                    xmax=[  x[1] for x in l_repr_ab_for_lollyplot ], 
                                    color='grey', alpha=.8)
                                    
        ax[1][1].plot([0,0],[min(my_range),max(my_range)], lw=0.5, c='red', linestyle = 'dotted')
        ax[1][2].plot([0,0],[min(my_range),max(my_range)], lw=0.5, c='red', linestyle = 'dotted')
        ax[1][1].scatter([ x[0] for x in l_repr_ab_for_lollyplot ], my_range, color='blue', alpha=1., label='repl')
        ax[1][2].scatter([ x[1] for x in l_repr_ab_for_lollyplot ], my_range, color='orange', alpha=1. , label='AB')
        ax[1][1].legend()
        ax[1][2].legend()
    else:
        print()
        print("For lollyplot choose one between I_want_log10,I_want_absolute,I_want_fractions")
        ax[1][1].text(0,0, "For lollyplot choose one between I_want_log10,I_want_absolute,I_want_fractions\n you chose nothing")
        print()
    if I_am_hacking_stuff:
        print("WARNING: this is just a little test, switch me off later!")
        ax[0][1].clear()
        ax[1][1].clear()
        ax[0][1].hlines(y=my_range, xmin=[ -x[0] for x in l_repr_ab_for_lollyplot ], 
                                    xmax=[  x[1] for x in l_repr_ab_for_lollyplot ], 
                                    color='grey', alpha=.8)
        ax[0][1].plot([0,0],[min(my_range),max(my_range)], lw=0.5, c='red', linestyle = 'dotted')
        ax[0][1].scatter([ -x[0] for x in l_repr_ab_for_lollyplot ], my_range, color='blue', alpha=1., label='repl')
        ax[0][1].scatter([  x[1] for x in l_repr_ab_for_lollyplot ], my_range, color='orange', alpha=1. , label='AB')
        ax[0][1].legend()
    
        sum_l_repr_ab_for_lollyplot_0 = sum(x[0] for x in l_repr_ab_for_lollyplot)
        sum_l_repr_ab_for_lollyplot_1 = sum(x[1] for x in l_repr_ab_for_lollyplot)
        print("sum_l_repr_ab_for_lollyplot_0,sum_l_repr_ab_for_lollyplot_1",sum_l_repr_ab_for_lollyplot_0,sum_l_repr_ab_for_lollyplot_1)
        print( [x[0] for x in l_repr_ab_for_lollyplot] )
        print( [x[1] for x in l_repr_ab_for_lollyplot] )
        # sys.exit(1)
        ax[0][1].clear()
        ax[1][1].clear()
        l_repr_ab_for_lollyplot = [ [x[0]/float(sum_l_repr_ab_for_lollyplot_0), x[1]/float(sum_l_repr_ab_for_lollyplot_1) ] for x in l_repr_ab_for_lollyplot ]
        # ax[1][1].hlines(y=my_range, xmin=[ -x[0] for x in l_repr_ab_for_lollyplot ], 
                                    # xmax=[  x[1] for x in l_repr_ab_for_lollyplot ], 
                                    # color='grey', alpha=.8)
        ax[1][1].plot([0,0],[min(my_range),max(my_range)], lw=0.5, c='red', linestyle = 'dotted')
        ax[1][1].scatter([ -x[0] for x in l_repr_ab_for_lollyplot ], my_range, color='blue', alpha=1., label='repl')
        ax[0][1].scatter([  x[1] for x in l_repr_ab_for_lollyplot ], my_range, color='orange', alpha=1. , label='AB')
        ax[1][1].legend()
        plt.show()
else:
    #Collect anyone that is alive now, sort by howmany ABs they made
    #then make cumulative plot
    l_ab_prod=[]
    l_ab_predicted=[]
    l_hist =[]
    for key in dg:
        #if dg[key]['death'] >= maxtime:
            #this guy is alive now
        ab_prod = dg[key]['A']
            #birthdate = dg[key]['birth']
            #l_ab_prod.append(ab_prod/float(maxtime-birthdate))
        l_ab_prod.append(ab_prod)
            #gnm = dg[key]['gnm']
            #nA = gnm.count('A') 
            #nF = gnm.count('F')
            #l_hist.append(nF)
            #l_ab_predicted.append( nA/(nA + 3.)*math.exp(-nF) )

    
    tot_people_here = len(l_ab_prod)
    print("tot people: ", tot_people_here)
    sum_ABs_here = sum(l_ab_prod)
    print("Sum ABs = ", sum_ABs_here)
    sorted_l_ab_prod = sorted(l_ab_prod, reverse=True)
    # print(len(l_ab_prod), sorted_l_ab_prod[0],sorted_l_ab_prod[-1] )
    from itertools import accumulate
    l_sum_sorted = list( accumulate(sorted_l_ab_prod) ) # cumulative sum A LOT faster than what I wrote!
    X = [x/float(tot_people_here) for x in range(tot_people_here)]
    Y = [y/float(sum_ABs_here) for y in l_sum_sorted]
    ax[1][1].plot( X , Y , label = 'cumulative antibiotic production' )
    
    # Attempt at log transforming the data
    # data is in [0,1], mapped to -log10(1-y)
    # result is nice and clear, but not really suited for anyting
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    #Create an inset in the lower right corner (loc=4) with borderpad=1, i.e.
    # 10 points padding (as 10pt is the default fontsize) to the parent axes
    axins = inset_axes(ax[1][1], width="100%", height="100%", bbox_to_anchor=(.6, .2, .4, .4), bbox_transform=ax[1][1].transAxes)

    lX=[]
    lY=[]
    for x,y in zip(X,Y):
        if y<0.999:
            try:
                lY.append(-math.log10(1.-y))
                lX.append(x)
            except:
                break
    axins.plot( lX , lY )

    for x,y in zip(X,Y):
        if y>=0.99:
            print("HEllo x,y", x,y)
            ax[1][1].scatter([x],[y],s=20, c='r')
            ax[1][1].plot([0,x,x],[y,y,0],ls='--', lw=0.5, c='r', label='99%')
            try:
                logy = -math.log10(1.-y)
                axins.scatter([x],[logy],s=10,c='r')
                axins.plot([0,x,x],[logy,logy,0],ls='--', lw=0.5, c='r')
            except ValueError:
                print("Math domain error, not plotting inset")
            break
    
    ax[1][1].legend()
    ax[1][1].set_xlim([0,1])
    ax[1][1].set_ylim([0,1])
    ax[1][1].set_xlabel('Cumulative fraction of population')
    ax[1][1].set_ylabel('Fraction of total antibiotic production')
    
    axins_yaxis = [0., 0.9,0.99,0.999]
    axins.set_yticks([-math.log10(1.-x) for x in axins_yaxis])
    axins.set_yticklabels([str(x) for x in axins_yaxis])
    # axins.set_y ([str(x) for x in axins_yaxis])
    
    #ax[1][1].plot( [x/float(tot_people_here) for x in range(tot_people_here)] ,[y/float(sum_ABs_here) for y in sorted_l_ab_prod])
    # ax[1][1].set_yscale('log')
    # mybins = np.linspace(0,1,21)
    # hist,bins  = np.histogram(l_ab_prod)
    # hist2,bins = np.histogram(l_ab_predicted, bins=mybins)
    
    # ax[1][1].plot(bins[:-1],hist, lw=1)
    # ax[1][1].plot(bins[:-1],hist2, lw=1, label='hist Ab predicted')
    #ax[1][1].plot( [ 0 , bins[-1] ] , [0,0], lw=0.1, c='red', linestyle = 'dotted' ) #a zero line to facilitate comparison
plt.show()
