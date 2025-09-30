#!/bin/bash
typeStim="b"
yposR=1.0
source_g=8000
cyt_conc=0.00005
dummy_spine_spacing=30.0
ns=2
ypos=0.25
coactive=1.0
nb=1
NstimL=10.0
freq=0.1
bg_type="b"
spacing=1.0
mixed=True
bgy=1
kcat_enz_tiam=20
for trial in 13
#for trial in 0 1 2 3 4 5 6 7 8 9
#for trial in 10 11 12 13 14 15 16 17 18 19
  do	
    for bfreq in 20.0
       do
          nohup python3 makeNewSpine.py --ns $ns --kagg -60.0 --mu0 -4.5 --spacing $spacing --nb $nb --freq $freq --kcat_enz_tiam $kcat_enz_tiam --yposR $yposR --ypos $ypos --bfreq $bfreq --s $source_g --typeStim $typeStim --spine_create_delay 5.0 --NstimL $NstimL --cyt_conc $cyt_conc  --dummy_spine_spacing $dummy_spine_spacing --coactive $coactive --trial $trial --mixed $mixed --bgy $bgy > ${kcat_enz_tiam}_${bgy}_${trial}.out&
       done
  done       
