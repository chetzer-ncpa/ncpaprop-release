set key top right samplen 2

set xlabel "Range [km]"
set ylabel "T.loss [kdB]" offset 0.25

set term pngcairo enh color font ",14" size 820,500
set out "../ePape_ex4_1d.png"

set title "1D Transmission Loss Magnitude; 0.5 Hz"

plot "code_outputs/ex4_tloss_1d_topo.pe" using 1:(20*log10(sqrt($3*$3+$4*$4))) lw 2 title "topo",\
	"code_outputs/ex4_tloss_1d_flat.pe" using 1:(20*log10(sqrt($3*$3+$4*$4))) lw 2 title "flat"

unset key

set xlabel "Range [km]"
set ylabel "Altitude [km]" offset 1
set cblabel "T.loss [dB re 1]" offset 1

#set palette defined (1 "white",2 "black")
#set palette defined (0 "white",1 "yellow",2 "red",3 "blue")
set palette defined ( 0 "light-blue", 1 "yellow", 2 "red" )

set title "2D Transmission Loss Magnitude topo; 0.5 Hz"
set size 0.95,1

set term pngcairo enh color font ",14" size 820,500
set out "../ePape_ex4_2d.png"

set cbrange [-130:-90]
plot [*:*] [-2:100] "code_outputs/ex4_tloss_2d_topo.pe" using 1:2:(20*log10(sqrt($3*$3+$4*$4))) palette,\
	"code_outputs/ex4_topography.pe" using (0.001*$2):(0.001*$3) with filledcurves x1 lc rgb '#784860'

