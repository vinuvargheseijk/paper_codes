#!/bin/bash
cytD=0.4
ns=1
delay=0.0
#for cyt_conc in 0.4e-3 0.42e-3 0.44e-3 0.46e-3 0.48e-3 0.5e-3 0.6e-3 0.7e-3 0.8e-3 1.0e-3 1.4e-3 1.8e-3 2.0e-3
#for cyt_conc in 0.6e-3 0.7e-3 0.8e-3 1.0e-3 1.4e-3 1.8e-3 2.0e-3
for cyt_conc in 0.6e-3 0.8e-3 1.2e-3 1.4e-3
  do
    for dummy_spine_spacing in 0.55
         do
          nohup python3 makeNewSpine.py --ns $ns --kagg -60.0 --mu0 -4.5 --cytD $cytD --dummy_spine_spacing $dummy_spine_spacing --cyt_conc $cyt_conc --spine_create_delay $delay > ${cyt_conc}_${dummy_spine_spacing}_${delay}.out&
         done
  done
