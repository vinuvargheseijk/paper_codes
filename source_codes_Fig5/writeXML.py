import os
import xml.etree.ElementTree as ET

def writeXML_td(time_d,fileName_td):
    if os.path.isfile(fileName_td):
       tree = ET.parse(fileName_td)
       root = tree.getroot()
       for i in range(len(time_d)):
          avec = ET.SubElement(root, 'time_d'+str(i)+str(fileName_td))
          avec.text = ''.join(str(j)+' ' for j in time_d[i])+'\n'
       tree.write(fileName_td)
    else:
       root = ET.Element('Data')
       for i in range(len(time_d)):
          avec = ET.SubElement(root, 'time_d'+str(i)+str(fileName_td))
          avec.text = ''.join(str(j)+' ' for j in time_d[i])+'\n'
       tree = ET.ElementTree(root)
       tree.write(fileName_td)
