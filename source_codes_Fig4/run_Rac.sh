#!/bin/bash

for ypos in 0.5e-6
  do	
    for src in 42000.0 52000.0
       do
          nohup python3 makeNewSpine.py --ns 1 --kagg -45.0 --mu0 -4.5 --spacing 0.75 --cytD 0.75 --nb 30 --phi_entire 6000.0 --freq 0.5 --freq_ratio 1.0 --ypos $ypos --bfreq 1.0 --bs 1.0 --s $src --CaTau 0.08 > n8000.out&
       done
  done       
