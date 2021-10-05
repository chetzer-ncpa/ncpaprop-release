#!/bin/bash 

gnuplot -e "dir='${NCPAPROP_DIR}/samples/' " profiles_ex2.gp 
gnuplot ePape_ex1.gp  
gnuplot ePape_ex2.gp  
gnuplot ePape_ex3.gp

