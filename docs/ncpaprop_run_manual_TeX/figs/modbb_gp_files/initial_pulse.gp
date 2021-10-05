
unset grid

set pointsize 2

set origin 0,0

set size 1.4,1
set term post enh eps color solid 22
set out "initial_pulse_model.eps"

set multiplot

set size 0.6,1
set origin 0,0
set title "Source Waveform"
set xlabel "Time [s]" 
set ylabel "Normalized Amplitude [ ]" offset 1
set key samplen 2
plot [0:6] [-1.5:1.5] "code_outputs/initial_pulse.dat" using 1:($2/3.23379) lw 4 title "s"

set size 0.8,1
set origin 0.6,0
set title "Source Fourier Transform"
set xlabel "Freq [Hz]"
set ylabel "Normalized Amplitude [1/Hz]" offset 1
set key samplen 2
plot "code_outputs/ft_initial_pulse.dat" using 1:($2/3.23379) lw 4 title "Re q",\
     "code_outputs/ft_initial_pulse.dat" using 1:($3/3.23379) lw 4 title "Im q",\
     "code_outputs/ft_initial_pulse.dat" using 1:(mag($2,$3)/3.23379) lw 4 title "|q|"

unset multiplot

! epstopdf initial_pulse_model.eps; rm initial_pulse_model.eps; mv initial_pulse_model.pdf ..

