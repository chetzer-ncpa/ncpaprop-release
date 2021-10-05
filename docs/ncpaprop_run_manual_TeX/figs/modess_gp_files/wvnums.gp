unset key
set style data lines

set style line 1 lc "black" lw 3
set style line 2 lc "dark-blue" lw 3
set style line 3 lc "dark-red" lw 3


set parametric

set xlabel "Soundspeed [m/s]
set ylabel "Altitude [km]" offset 1

set arrow from 208,94.5 to 453,94.5 heads ls 2
set arrow from 314,9.2 to 453,9.2 heads ls 3
set title "Relevant Phase Speeds"

set term pngcairo enh color font ",14" size 600,500
set out "../wvnums_modess.png"

plot [0:20] [150:500] [0:160] "code_outputs/toymodel_60_ms_jet.dat" using (20*sqrt($5)+$2):1 ls 1,\
     206.4,t+84.5 ls 2,454,t+84.5 ls 2,\
     314,t ls 3 ,454,t ls 3

unset arrow

