unset key 

set pm3d map
set palette defined (0 "white",1 "yellow",2 "red",3 "blue")

unset grid
set angles radians

c(x,y)=x*cos(pi*y/180.0)
s(x,y)=x*sin(pi*y/180.0)

set xlabel "Easterly [km]" 
set ylabel "Northerly [km]"
set xtics rotate
set cblabel "Tloss [dB]" offset 1

set xrange [-1000:1000]
set yrange [-1000:1000]
set cbrange [-130:-90]

dir="code_outputs/"

set term pngcairo enh color font ",14" size 1200,600
set out "../ePape_ex3.png"

set multiplot
set size 0.45,0.9
set size square

set origin 0.02,0.1
set title "Transmission Loss Magnitude; 0.1 Hz lossy"
splot dir."ex3_tloss_multiprop.pe" using (s($1,$2)):(c($1,$2)):(20*log10(sqrt($3*$3+$4*$4)))

set origin 0.52,0.1
set title "Transmission Loss Magnitude; 0.1 Hz lossless"
splot dir."ex3_tloss_multiprop_lossless.pe" using (s($1,$2)):(c($1,$2)):(20*log10(sqrt($3*$3+$4*$4)))

unset multiplot

