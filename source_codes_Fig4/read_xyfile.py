
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import sys
import xml.etree.ElementTree as ET
from matplotlib import animation
import math
import totalEnergy as TE

fig, (ax1, ax2) = plt.subplots(2, 1)
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

ax1.set_xlim([1, 10])
ax1.set_ylim([0, 0.15])
ax1.set_xticks([])
ax1.set_yticks([])

ax2.set_xlim([1, 10])
ax2.set_ylim([0, 5])
ax2.set_xticks([])
ax2.set_yticks([])

ax1.set_xlabel("$\mu m$")
ax1.set_ylabel("Height $\mu m$")

ax2.set_xlabel("$\mu m$")
ax2.set_ylabel("Cytosol")

line1, = ax1.plot([], [], lw=2)
line2, = ax2.plot([], [], lw=2)
line1.set_color('blue')
line_in, = ax1.plot([], [], lw=2)
line_in.set_color('black')
texts_in=ax1.text(.2,.8, "{}".format(0.1), transform=ax1.transAxes, color = 'k')
texts=ax1.text(.6,.8, "{}".format(0.1), transform=ax1.transAxes, color = 'b')
def plotXML(filename,tree):
   an_y=[]
   max_an_y=[]
   sum_an=[]
   ind_hs=[]
   ind_fs=[]
   for k in tree:
    m = 0   
    break_flag = 0
    while (break_flag ==0):
       yValues = []
       try:
          yaxis = k.find("time_d"+str(m)+filename)
          yValues = [float(j) for j in yaxis.text.split()]
       except:
          break_flag = 1
       if break_flag == 0:   
         an_y.append(yValues)
         max_an_y.append(max(yValues))
         sum_an.append(sum(yValues))
       m = m + 20

   return sum_an,max_an_y,an_y

def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line_in.set_data([], [])
    return line1, line_in, line2
   

def animate(i):
    #x = np.asarray(an_yX[i])
    x = np.asarray(an_yX[i])
    y = np.asarray(an_yY[i])
    c = np.asarray(an_yC[i])
    #texts.set_text('t= '+str(i * 0.25))
    texts.set_text('Evolving')
    texts_in.set_text('Initial')
    line1.set_data(x, y)
    line2.set_data(x_vec, c)
    plt.tight_layout()
    
    #line_in.set_data(an_yX[0],an_yY[0])
    return line1,line_in, line2


num_spines = 3
strength = 8000
spacing = 2.0
scaling = 0.4
directory = "./"
   

xmlfilename1 = "X.xml"
xmlfilename2 = "Y.xml"
xmlfilename3 = "cytosol.xml"

xmlfiles = [xmlfilename1,xmlfilename2, xmlfilename3]
for root, dirs, files in os.walk(directory):
    for xmlfile in range(len(xmlfiles)): 
        for file in files:  
            if file == xmlfiles[xmlfile]:
               print("File under processing: ", xmlfiles[xmlfile]) 
               tree1 = ET.parse(directory + "/" + xmlfiles[xmlfile])
               tree = [tree1]
               if xmlfiles[xmlfile] == xmlfilename1:
                 sum_anX,max_an_yX,an_yX=plotXML(xmlfiles[xmlfile],tree)
               if xmlfiles[xmlfile] == xmlfilename2:
                 sum_anY,max_an_yY,an_yY=plotXML(xmlfiles[xmlfile],tree)
               if xmlfiles[xmlfile] == xmlfilename3:
                 sum_anC,max_an_yC,an_yC=plotXML(xmlfiles[xmlfile],tree)
                 num_voxels = len(an_yC[-1])
                 print("Length: ", TE.Length)
                 x_vec = np.linspace(0, TE.Length, num_voxels)


an_yX_scaled = []
for i in range(len(an_yX)):
  an_yX_scaled.append(np.asarray(an_yX[i]) * 2.5)

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(an_yY), interval=25, blit=True)

anim.save('shapeN3.mp4', fps=25, extra_args=['-vcodec', 'libx264'])

plt.show()
