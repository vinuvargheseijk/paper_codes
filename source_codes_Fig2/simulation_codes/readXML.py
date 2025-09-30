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
       m = m + 1

   return sum_an,max_an_y,an_y

