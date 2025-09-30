#!/bin/bash

mu0=-4.5
spacing=2.0
dt=0.24
cytD=0.4

for agg in -45.0 -50.0 -55.0 -60.0 -65.0
  do	
    for phi_entire in 2000.0 3000.0 4000.0 5000.0 6000.0 7000.0 8000.0
       do
          nohup python3 makeNewSpine.py --phi_entire $phi_entire --ns 1.0 --kagg $agg --mu0 $mu0 --spacing $spacing --dt $dt --cytD $cytD > n_${phi_entire}_${spacing}_${ns}_${k_agg}.out &
       done
  done       
