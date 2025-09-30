import numpy as np
import os
import pandas as pd

max_vars = 10
def extract_numbers(file):
    if len(file) > 10:
        extension = file[-7] + file[-6] + file[-5] + file[-4] + file[-3] + file[-2] + file[-1]
        numbers = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", ".", "-"]
        var_count = "v"
        var_name = "n"
        vars_values = []
        vars_names = []
        for var in range(max_vars):
            globals()[var_count] = ""
            globals()[var_name] = ""
            first_val_occurance = 0
            for i in file:
              if i in numbers:
                 globals()[var_count] = globals()[var_count] + i
                 first_val_occurance = 1
              if i not in numbers and first_val_occurance == 1:
                 vars_values.append(globals()[var_count])
                 vars_names.append(globals()[var_name])
                 break
              if i not in numbers and first_val_occurance == 0:
                if i != "_": 
                  globals()[var_name] = globals()[var_name] + i
               
              file = file.replace(i, "", 1)    
    else:
        extension = "File not needed"
        vars_values = []
        vars_names = []
    return extension, vars_values, vars_names       


filetype = "csv"

files = os.listdir("./")

for i in range(max_vars):
    globals()["var" + str(i)] = []

max_size_list = 0
for f in files:
       extension, vars_values, vars_names = extract_numbers(f)
       if len(vars_values) > max_size_list:
           max_size_list = len(vars_values)
       if extension in "phi.csv":
          for vc in range(len(vars_values)):
              try:  
                 globals()["var" + str(vc)].append(float(vars_values[vc]))
              except:
                 print("Not a variable")


df = pd.DataFrame()
for i in range(max_size_list - 1):
    print(globals()["var" + str(i)])
    df[vars_names[i]] = globals()["var" + str(i)]

print(df)    




