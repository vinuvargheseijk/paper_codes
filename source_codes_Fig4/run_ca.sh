#!/bin/bash

for ypos in 0.5e-6 1.0e-6 1.5e-6 2.0e-6
  do	
    for catau in 0.02 0.04 0.06 0.08
       do
          nohup python3 makeNewSpine.py --ns 1 --kagg -45.0 --mu0 -4.5 --spacing 0.75 --cytD 0.75 --nb 1 --phi_entire 6000.0 --freq 1.0 --freq_ratio 1.0 --ypos $ypos --bfreq 1.0 --bs 1.0 --s 16000.0 --CaTau $catau > n8000.out&
       done
  done       
