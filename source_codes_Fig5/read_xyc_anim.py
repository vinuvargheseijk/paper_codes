import tag
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import sys
import xml.etree.ElementTree as ET
from matplotlib import animation
import math
import totalEnergy as TE
import readXML
import pandas as pd

fig = plt.figure(figsize = (8,6))
ax1 = plt.axes( xlim = (0, TE.Length), ylim = (0, 0.15) )
ax2 = ax1.twinx()
file_locations = pd.read_csv(tag.tag_string + "MaxLocations.csv")
total_spines = len(file_locations)

#ax1.set_xlim([3.5, 6.5])
ax1.set_ylim([0, 0.5])

ax2.set_ylim([0, 10.0])

ax1.set_xlabel("$\mu m$")
ax1.set_ylabel("Height $\mu m$")

ax2.set_ylabel("Cytosol")

line1, = ax1.plot([], [], lw=2)
line2, = ax2.plot([], [], lw=2)
line3, = ax2.plot([], [], lw=2)
line4, = ax2.plot([], [], lw=2)
line1.set_color('blue')
line2.set_color('red')
line3.set_color('green')
line4.set_color('magenta')
line_in, = ax1.plot([], [], lw=2)
line_in.set_color('magenta')
texts_in=ax1.text(.2,.8, "{}".format(0.1), transform=ax1.transAxes, color = 'k')
texts=ax1.text(.6,.8, "{}".format(0.1), transform=ax1.transAxes, color = 'b')
texts_cyt=ax2.text(.4,.8, "{}".format(0.1), transform=ax1.transAxes, color = 'g')
texts_ca=ax2.text(.6,.4, "{}".format(0.1), transform=ax1.transAxes, color = 'r')
texts_mem=ax2.text(.8,.8, "{}".format(0.1), transform=ax1.transAxes, color = 'magenta')

def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    line_in.set_data([], [])
    return line1, line_in, line2, line3
   
def save_total_conc(var1, var2):
    for spines in range(1,total_spines):
        spine_index_to_plot = spines
        start_of_the_spine = file_locations.iloc[spine_index_to_plot]["dome_start"] - 0.2e-6
        start_of_the_spine = int(start_of_the_spine / TE.diffL) 
        end_of_the_spine = file_locations.iloc[spine_index_to_plot]["dome_end"] + 0.2e-6
        end_of_the_spine = int(end_of_the_spine / TE.diffL) 
        globals()["total_conc_save" + str(spines)] = []
        globals()["marker" + str(spines)] = []
        for i in range(len(var1)):
           sum_mem = sum(var2[i][start_of_the_spine : end_of_the_spine])
           sum_cyt = sum(var1[i])
           total = sum_mem + sum_cyt
           total_conc = 1e3 * total / (TE.V_cyt * TE.Na)
           globals()["total_conc_save" + str(spines)].append(total_conc)
           if sum_mem > 1e-2:
              globals()["marker" + str(spines)].append(1)
           else:
              globals()["marker" + str(spines)].append(0)
    all_spine_conc = []
    all_spine_marker = []
    for spines in range(1, total_spines):
        all_spine_conc.append(globals()["total_conc_save" + str(spines)])
        all_spine_marker.append(globals()["marker" + str(spines)])
    
    return all_spine_conc, all_spine_marker


def animate(i):
    time = chemTimes[i]
    if len(shapeTimes) > 0:
      if time > shapeTimes[0]:
          x = np.asarray(an_yX[0])
          y = np.asarray(an_yY[0])
          line1.set_data(x, y)
          an_yX.pop(0)
          an_yY.pop(0)
          shapeTimes.pop(0)
    else:  
      line1.set_data(0, 0)
    ca = np.asarray(an_yCa[i])
    cyt = np.asarray(an_yCyt[i])
    mem = np.asarray(an_yMem[i])
    ca_conc = 1e3 * ca / (TE.Na * TE.V_voxel)
    cyt_conc = 1e3 * cyt / (TE.Na * TE.V_voxel)
    mem_conc = 1e3 * mem / (TE.Na * TE.V_voxel)
    total = np.asarray(cyt) + np.asarray(mem)
    total_conc = 1e3 * total / (TE.Na * TE.V_voxel)
    texts_in.set_text('t = '+ str(round(time,1)))
    texts.set_text('Evolving shape')
    #texts_in.set_text('Initial')
    texts_ca.set_text('Calcium')
    texts_cyt.set_text('Cytosol')
    texts_mem.set_text('Total(Mem + Cyt)')
    line2.set_data(x_vec, ca_conc)
    line3.set_data(x_vec, cyt_conc)
    line4.set_data(x_vec, total_conc)
    plt.tight_layout()
    
    #line_in.set_data(an_yX[0],an_yY[0])
    return line1,line_in, line2, line3, line4

directory = "./"
tag_string = tag.tag_string



xmlfilename1 = tag_string+"X.xml"
xmlfilename2 = tag_string+"Y.xml"
xmlfilename3 = tag_string+"Ca.xml"
xmlfilename4 = tag_string+"cytosol.xml"
xmlfilename5 = tag_string+"membrane.xml"

shapeTimes = list(pd.read_csv(tag_string + "shapeTimes.csv")["0"])
chemTimes = list(pd.read_csv(tag_string + "chemTimes.csv")["time"])

xmlfiles = [xmlfilename1,xmlfilename2, xmlfilename3, xmlfilename4, xmlfilename5]
for root, dirs, files in os.walk(directory):
    for xmlfile in range(len(xmlfiles)): 
        for file in files:  
            if file == xmlfiles[xmlfile]:
               print("File under processing: ", xmlfiles[xmlfile]) 
               tree1 = ET.parse(directory + "/" + xmlfiles[xmlfile])
               tree = [tree1]
               if xmlfiles[xmlfile] == xmlfilename1:
                 sum_anX,max_an_yX,an_yX=readXML.plotXML(xmlfiles[xmlfile],tree)
                 print(len(an_yX))
               if xmlfiles[xmlfile] == xmlfilename2:
                 sum_anY,max_an_yY,an_yY=readXML.plotXML(xmlfiles[xmlfile],tree)
                 print(len(an_yY))
               if xmlfiles[xmlfile] == xmlfilename3:
                 sum_anCa,max_an_yCa,an_yCa=readXML.plotXML(xmlfiles[xmlfile],tree)
                 print(len(an_yCa))
                 num_voxels = len(an_yCa[-1])
                 print("Length: ", TE.Length)
                 x_vec = np.linspace(0, TE.Length, num_voxels)
               if xmlfiles[xmlfile] == xmlfilename4:
                 sum_anCyt,max_an_yCyt,an_yCyt=readXML.plotXML(xmlfiles[xmlfile],tree)
                 print(len(an_yCyt))
               if xmlfiles[xmlfile] == xmlfilename5:
                 sum_anMem,max_an_yMem,an_yMem=readXML.plotXML(xmlfiles[xmlfile],tree)
                 print(len(an_yMem))


anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(an_yCa), interval=1, blit=True)

all_spine_conc, all_spine_marker = save_total_conc(an_yCyt, an_yMem)
print("Total conc: ", len(all_spine_conc))
df_temporal_conc = pd.DataFrame()
df_marker = pd.DataFrame()
for i in range(1, total_spines):
   print(i) 
   df_temporal_conc["conc" + str(i)] = all_spine_conc[i - 1]
   df_marker["conc" + str(i)] = all_spine_marker[i - 1]
df_temporal_conc.to_csv(tag.tag_string + "saved_conc.csv")
df_marker.to_csv(tag.tag_string + "saved_marker.csv")

anim.save(tag_string + 'shapeN3.mp4', fps=25, extra_args=['-vcodec', 'libx264'])

