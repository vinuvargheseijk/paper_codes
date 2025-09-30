import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import tag
import os
import totalEnergy as TE
import xml.etree.ElementTree as ET
import readXML
import sys


plot_interval = 10
dt = 0.02


def vol_dist():
    taper = tag.taper
    Length = tag.Length * 1e-6
    dendDia = TE.dendDia
    Na = TE.Na
    num_voxles = Length / TE.diffL
    comptLen = 1e-6
    comptLen_voxel = comptLen / TE.diffL
    vol_list = []
    conc_fact = []
    segn = 0
    for num in range(num_voxels):
        if num % comptLen_voxel == 0:
           dia_seg = dendDia - segn * taper * dendDia
           vol_seg = np.pi * (0.5 * dia_seg)**2 * TE.diffL
           vol_list.append(vol_seg)
           conc_fact.append(1 / (vol_seg * Na))
           segn = segn + 1
        else:   
           print("Dia: ", dia_seg)
           print("Volume: ", vol_seg)
           vol_list.append(vol_seg)
           conc_fact.append(1 / (vol_seg * Na))
    print("Num voxels, num segments: ", num_voxels, len(vol_list))       
    return vol_list, conc_fact        


def plot_radar(filename, axis):
  df = pd.read_csv("./energies1.csv")
  categories = ["Ae", "Ee", "Me", "Ae"]
  c = list(df.iloc[0][1:])
  print(c)
  print(len(c))
  c.append(c[0])
  angles = np.linspace(0, 2 * np.pi, len(c))
  print(angles)
  axis.plot(angles, c)
  axis.fill(angles, c)
  lines, labels = axis.set_thetagrids(np.degrees(angles), labels=categories)
  return lines, labels

def find_bigger_df(df1, df2):
    if len(df1) > len(df2):
        big = df1
    else:
        big = df2
    return big 


#fig, (ax1, ax2) = plt.subplots(2, subplot_kw={'projection': 'polar'})
fig = plt.figure(figsize = (10,8))
ax = fig.add_subplot(221, projection = 'polar')
ax1 = fig.add_subplot(222, projection = 'polar')
ax2 = fig.add_subplot(223)
ax.set_ylim(-10000, 10000)
ax1.set_ylim(-10000, 10000)
#ax2.set_ylim(0, 0.003)
ax2.set_xlim(0, tag.Length)
ax3 = ax2.twinx()
ax2.set_ylim(0,0.05)
ax3.set_ylim(0, 5e-3)
texts_in=ax3.text(.2,-0.6, "{}".format(0.1), transform=ax1.transAxes, color = 'k')
texts=ax2.text(.1,-0.7, "{}".format(0.1), transform=ax1.transAxes, color = 'b')
texts_cyt=ax3.text(.1,-0.8, "{}".format(0.1), transform=ax1.transAxes, color = 'g')
texts_cytC=ax3.text(.1,-0.9, "{}".format(0.1), transform=ax1.transAxes, color = 'r')
ax3.set_xlim(0, tag.Length)
ax2.set_xlabel("$\mu m$")
ax2.set_ylabel("Height $\mu m$")
"""
plot_radar("./energies1.csv", ax1)
plot_radar("./energies2.csv", ax2)
"""
directory = "./"
tag_string = tag.tag_string

xmlfilename1 = tag_string+"X.xml"
xmlfilename2 = tag_string+"Y.xml"
xmlfilename3 = tag_string+"marker.xml"
xmlfilename4 = tag_string+"cytosolConc.xml"

filename1 = "./energies1.csv"
filename2 = "./energies2.csv"
l, = ax.plot([],[])
l1, = ax1.plot([],[])
l2, = ax2.plot([],[])
l3, = ax3.plot([],[])
l4, = ax2.plot([],[])
l3.set_color("green")
l4.set_color("red")
df1 = pd.read_csv(filename1)
df2 = pd.read_csv(filename2)
categories = ["Aggregation", "Entropy", "Mismatch", "Bending", "Aggregation"]
big_df = find_bigger_df(df1, df2)




def update(i):
  x = np.asarray(an_yX[i])
  y = np.asarray(an_yY[i])
  cyt = np.asarray(an_yC[i])
  m = np.asarray(an_yM[i])
  ax.clear()  
  ax1.clear()
  global df1  
  global df2 
  c = list(df1.iloc[i][1:])
  c1 = list(df2.iloc[i][1:])
  c.append(c[0])
  c1.append(c1[0])
  texts.set_text("Shape")
  texts_cyt.set_text("Cyt Conc")
  texts_cytC.set_text("Marker")
  angles = np.linspace(0, 2 * np.pi, len(c))
  l.set_data(angles, c)
  l1.set_data(angles, c1)
  l2.set_data(x, y)
  vol_list, conc_fact = vol_dist()
  cyt_conc = np.asarray(cyt) * conc_fact
  l3.set_data(x_vec,m)
  l4.set_data(x_vec, cyt)
  texts_in.set_text('t= '+str(i * dt * plot_interval))
  ax.set_thetagrids(np.degrees(angles), labels=categories)
  ax1.set_thetagrids(np.degrees(angles), labels=categories)
  if c[0] != 0:
    ax.fill(angles, c)
  if c1[0] != 0:  
    ax1.fill(angles, c1)
  return l,
    
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
                 x_vec = np.linspace(0, TE.Length, num_voxels)
                 print("Length: ", TE.Length)
               if xmlfiles[xmlfile] == xmlfilename4:
                 sum_anM,max_an_yM,an_yM=readXML.plotXML(xmlfiles[xmlfile],tree)
                 print(len(an_yM))

ani = animation.FuncAnimation(fig, update, frames = len(an_yM), interval = 1, blit = True)
#ani = animation.FuncAnimation(fig, update, frames = 100, interval = 1, blit = True)
ani.save('shapeN3.mp4', fps=25, extra_args=['-vcodec', 'libx264'])
#plt.show()



