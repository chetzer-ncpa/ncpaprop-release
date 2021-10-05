unset key 

mag(x,y)=sqrt(x**2+y**2)
TL(x)=20.0*log10(x)

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

set term pngcairo enh color font ",14" size 820,500
set size 0.95,1

set out "../wmod_ex2.png"
splot "code_outputs/wtloss_2d.nm" using 1:2:(TL(mag($3,$4)))


