set key samplen 2

set style line 1 linecolor "black" lw 2 

dir="code_outputs/"

set xlabel "Range [km]"
set ylabel "TL [dB]"

set title "1D Transmission Loss Magnitude; 0.5 Hz"

set term pngcairo enh color font ",14" size 820,500

set out "../ePape_ex2_1d.png"
plot [0:1200] [*:*] dir."ex2_tloss_1d.pe" using 2:(20*log10(sqrt($3*$3+$4*$4))) ls 1 title "" with lines

unset key 

#set palette defined (0 "white",1 "yellow",2 "red",3 "blue")
set palette defined ( 0 "light-blue", 1 "yellow", 2 "red" )


set xtics 100
set ytics 50
set cbtics 10
set xrange [0:1200]
set xtics 200
set yrange [0:120]
set ytics 30
set cbrange [-130:-90]
set xlabel "Range [km]"
set ylabel "Altitude [km]" offset 0.5
set cblabel "Tloss [dB]" offset 1

set title "2D Transmission Loss Magnitude; 0.5 Hz"
set size 0.95,1

set out "../ePape_ex2_2d.png"
plot [*:*] [-2:100] dir."ex2_tloss_2d_flat.pe" using 1:2:(20*log10(sqrt($3*$3+$4*$4))) palette,\
	0 with filledcurves x1 lc rgb '#784860'


