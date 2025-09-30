#!/bin/bash
ypos=0.5e-6
for nb in 30 45
  do	
    for src in 32000
       do
          nohup python3 makeNewSpine.py --ns 1 --kagg -45.0 --mu0 -4.5 --spacing 0.75 --cytD 0.75 --nb $nb --phi_entire 6000.0 --freq 0.5 --freq_ratio 1.0 --ypos $ypos --bfreq 1.0 --bs 1.0 --s $src --CaTau 0.08 > n8000.out&
       done
  done       
