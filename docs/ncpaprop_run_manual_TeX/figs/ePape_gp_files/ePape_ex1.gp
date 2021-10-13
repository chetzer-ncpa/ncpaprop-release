set key samplen 2

dir="code_outputs/"
TL(x)=20*log10(x)

set xlabel "Range [km]"
set ylabel "TL [dB]"

set title "1D Transmission Loss Magnitude; 0.1 Hz"

set term pngcairo enh color font ",14" size 820,500

set out "../ePape_ex1_1d.png"
plot dir."ex1_tloss_1d.pe" using 1:(TL(mag($3,$4))) lt 6 lw 3 title "" with lines

unset key 

set pm3d map
set palette defined (0 "white",1 "yellow",2 "red",3 "blue")


set xtics 100
set ytics 50
set cbtics 10
set xrange [0:1000]
set yrange [0:150]
set cbrange [-130:-90]
set xlabel "Range [km]"
set ylabel "Altitude [km]"
set cblabel "Tloss [dB]" offset 1

set title "2D Transmission Loss Magnitude; 0.1 Hz"
set size 0.95,1

set out "../ePape_ex1_2d.png"
splot dir."ex1_tloss_2d.pe" using 1:2:(20*log10(sqrt($3*$3+$4*$4)))

