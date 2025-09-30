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
import seaborn as sns
sns.set_context("poster")
fig = plt.figure(figsize = (12,10))
#ax1 = plt.axes( xlim = (0, TE.Length), ylim = (0, 0.15) )
ax1 = fig.add_subplot(221)
ax1.set_xlim(0, TE.Length)
ax1.set_ylim(0, 0.15)
ax2 = ax1.twinx()
ax3 = fig.add_subplot(223, projection = 'polar')
ax4 = fig.add_subplot(224, projection = 'polar')
ax3.set_ylim(-15000, 6000)
ax4.set_ylim(-15000, 6000)
energies1 = pd.read_csv("./" + tag.tag_string + "energies1.csv")
energies2 = pd.read_csv("./" + tag.tag_string + "energies2.csv")
AE_en1 = list(energies1["Ae1"])
AE_en2 = list(energies2["Ae2"])

EE_en1 = list(energies1["Ee1"]) 
EE_en2 = list(energies2["Ee2"])


ME_en1 = list(energies1["Me1"])
ME_en2 = list(energies2["Me2"])

BE_en1 = list(energies1["Be1"]) 
BE_en2 = list(energies2["Be2"])

CE_en1 = list(energies1["Ce1"]) 
CE_en2 = list(energies2["Ce2"])

#ax1.set_xlim([3.5, 6.5])
ax1.set_ylim([0, 1.0])

ax2.set_ylim([0, 10.0])

ax1.set_xlabel("$\mu m$")
ax1.set_ylabel("Height $\mu m$")

ax2.set_ylabel("Cytosol")

categories = ["Aggregation", "Entropy", "Mismatch", "ChemP", "Aggregation"]

line1, = ax1.plot([], [], lw=2)
line2, = ax2.plot([], [], lw=2)
line3, = ax2.plot([], [], lw=2)
line4, = ax2.plot([], [], lw=2)
lineE1, = ax3.plot([], [], lw=2)
lineE2, = ax4.plot([], [], lw=2)
line1.set_color('blue')
line2.set_color('red')
line3.set_color('green')
line4.set_color('magenta')
line_in, = ax1.plot([], [], lw=2)
line_in.set_color('magenta')
texts_in=ax1.text(.2,.8, "{}".format(0.1), transform=ax1.transAxes, color = 'k')
texts=ax1.text(.6,.8, "{}".format(0.1), transform=ax1.transAxes, color = 'b')
texts_m=ax2.text(.4,.8, "{}".format(0.1), transform=ax1.transAxes, color = 'g')
texts_c=ax2.text(.6,.4, "{}".format(0.1), transform=ax1.transAxes, color = 'r')
texts_mem=ax2.text(.6,.2, "{}".format(0.1), transform=ax1.transAxes, color = 'magenta')

def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    lineE1.set_data([], [])
    lineE2.set_data([], [])
    line_in.set_data([], [])
    return line1, line_in, line2, line3, lineE1, lineE2
   


def animate(i):
    time = chemTimes[i]
    if time > shapeTimes[0]:
          x = np.asarray(an_yX[0])
          y = np.asarray(an_yY[0])
          line1.set_data(x, y)
          an_yX.pop(0)
          an_yY.pop(0)
          shapeTimes.pop(0)
          """
          e1 = list(energies1.iloc[0][1:])      
          e2 = list(energies2.iloc[0][1:])   
          e1.append(e1[0])
          e2.append(e2[0])
          """
          e1 = [AE_en1[0], EE_en1[0], ME_en1[0], CE_en1[0], AE_en1[0]]
          e2 = [AE_en2[0], EE_en2[0], ME_en2[0], CE_en2[0], AE_en2[0]]
          AE_en1.pop(0)
          EE_en1.pop(0)
          ME_en1.pop(0)
          BE_en1.pop(0)
          CE_en1.pop(0)

          AE_en2.pop(0)
          EE_en2.pop(0)
          ME_en2.pop(0)
          BE_en2.pop(0)
          CE_en2.pop(0)

          angles = np.linspace(0, 2 * np.pi, len(e1))
          lineE1.set_data(angles, e1)
          lineE2.set_data(angles, e2)
          ax3.set_thetagrids(np.degrees(angles), labels=categories)
          ax4.set_thetagrids(np.degrees(angles), labels=categories)
    c = np.asarray(an_yC[i])
    m = np.asarray(an_yM[i])
    mem = np.asarray(an_yMem[i])
    c = 1e3 * c / (TE.Na * TE.V_voxel)
    m = 1e3 * m / (TE.Na * TE.V_voxel)
    mem = 1e3 * mem / (TE.Na * TE.V_voxel)
    texts_in.set_text('t = '+ str(round(time,1)))
    texts.set_text('Evolving shape')
    #texts_in.set_text('Initial')
    texts_c.set_text('Calcium')
    texts_m.set_text('Cytosol')
    texts_mem.set_text('Membrane')
    line2.set_data(x_vec, c)
    line3.set_data(x_vec, m)
    line4.set_data(x_vec, mem)
    plt.tight_layout()
    
    #line_in.set_data(an_yX[0],an_yY[0])
    return line1,line_in, line2, line3, line4, lineE1, lineE2

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
                 sum_anC,max_an_yC,an_yC=readXML.plotXML(xmlfiles[xmlfile],tree)
                 print(len(an_yC))
                 num_voxels = len(an_yC[-1])
                 print("Length: ", TE.Length)
                 x_vec = np.linspace(0, TE.Length, num_voxels)
               if xmlfiles[xmlfile] == xmlfilename4:
                 sum_anM,max_an_yM,an_yM=readXML.plotXML(xmlfiles[xmlfile],tree)
                 print(len(an_yM))
               if xmlfiles[xmlfile] == xmlfilename5:
                 sum_anMem,max_an_yMem,an_yMem=readXML.plotXML(xmlfiles[xmlfile],tree)
                 print(len(an_yM))



anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(an_yC), interval=1, blit=True)
#anim = animation.FuncAnimation(fig, animate, init_func=init, frames=100, interval=1, blit=True)


anim.save(tag_string + 'shapeN3.mp4', fps=25, extra_args=['-vcodec', 'libx264'])

