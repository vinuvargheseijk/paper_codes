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

fig = plt.figure(figsize = (8,6))
ax1 = plt.axes( xlim = (0, TE.Length), ylim = (0, 0.15) )
ax2 = ax1.twinx()

#ax1.set_xlim([3.5, 6.5])
ax1.set_ylim([0, 1.0])

ax2.set_ylim([0, 30])

ax1.set_xlabel("$\mu m$")
ax1.set_ylabel("Height $\mu m$")

ax2.set_ylabel("Cytosol")

line1, = ax1.plot([], [], lw=2)
line2, = ax2.plot([], [], lw=2)
line3, = ax2.plot([], [], lw=2)
line1.set_color('blue')
line2.set_color('red')
line3.set_color('green')
line_in, = ax1.plot([], [], lw=2)
line_in.set_color('black')
texts_in=ax1.text(.2,.8, "{}".format(0.1), transform=ax1.transAxes, color = 'k')
texts=ax1.text(.6,.8, "{}".format(0.1), transform=ax1.transAxes, color = 'b')
texts_m=ax1.text(.4,.8, "{}".format(0.1), transform=ax1.transAxes, color = 'g')
texts_c=ax1.text(.6,.4, "{}".format(0.1), transform=ax1.transAxes, color = 'r')

def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    line_in.set_data([], [])
    return line1, line_in, line2, line3
   

def animate(i):
    x = np.asarray(an_yX[i])
    y = np.asarray(an_yY[i])
    c = np.asarray(an_yC[i])
    m = np.asarray(an_yM[i])
    texts_in.set_text('t= '+str(i * 0.025 * 10))
    texts.set_text('Evolving shape')
    #texts_in.set_text('Initial')
    texts_c.set_text('membrane')
    texts_m.set_text('marker')
    line1.set_data(x, y)
    line2.set_data(x_vec, c)
    line3.set_data(x_vec, m)
    plt.tight_layout()
    
    #line_in.set_data(an_yX[0],an_yY[0])
    return line1,line_in, line2, line3

directory = "./"
tag_string = tag.tag_string



xmlfilename1 = tag_string+"X.xml"
xmlfilename2 = tag_string+"Y.xml"
xmlfilename3 = tag_string+"membrane.xml"
xmlfilename4 = tag_string+"marker.xml"

xmlfiles = [xmlfilename1,xmlfilename2, xmlfilename3, xmlfilename4]
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



an_yX_scaled = []
for i in range(len(an_yX)):
  an_yX_scaled.append(np.asarray(an_yX[i]) * 2.5)

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(an_yC), interval=1, blit=True)

anim.save('shapeN3.mp4', fps=25, extra_args=['-vcodec', 'libx264'])

plt.show()

"""
plt.figure(101)
plt.plot(sum_anM, label = 'membrane')
plt.plot(sum_anC, label = 'cytosol')

plt.figure(102)
plt.plot(np.asarray(sum_anM) + np.asarray(sum_anC), label = 'sum')
plt.legend()
plt.show()
"""
