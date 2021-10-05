
unset grid
set origin 0,0

scale=1.68554

set size 1.4,1
set term post enh eps color solid 22
set out "example_source_model.eps"

set multiplot

set size 0.6,1
set origin 0,0
set title "Source Waveform"
set xlabel "Time [s]" 
set ylabel "Normalized Amplitude [ ]" offset 1
set key samplen 2
plot [0:25] [-1.5:1.5] "code_outputs/source_waveform.dat" using 1:($2/scale) lw 4 title "s"

set size 0.8,1
set origin 0.6,0
set title "Source Fourier Transform"
set xlabel "Freq [Hz]"
set ylabel "Normalized Amplitude [1/Hz]" offset 1
set key samplen 2
plot "code_outputs/source_spectrum.dat" using 1:($2/scale) lw 4 title "Re q",\
     "code_outputs/source_spectrum.dat" using 1:($3/scale) lw 4 title "Im q",\
     "code_outputs/source_spectrum.dat" using 1:(mag($2,$3)/scale) lw 4 title "|q|"

unset multiplot

! epstopdf example_source_model.eps; rm example_source_model.eps
! mv example_source_model.pdf ../modbb_ex_source_model.pdf

unset key
set origin 0,0
set size 1,0.8
set term post enh eps color solid 22
set out "example_prop_pulse.eps"

set title "Propagated Waveform at 270 km"
set xlabel "Travel Time [s]" 
set ylabel "Normalized Amplitude 10^{-7} [ ]" offset 1
plot [880:960] "code_outputs/mywavf.dat" using 1:(1.0e7*$2/scale) lt 3 lw 4 

! epstopdf example_prop_pulse.eps; rm example_prop_pulse.eps
! mv example_prop_pulse.pdf ../modbb_ex_prop_pulse.pdf
! cp ../modbb_ex_prop_pulse.pdf ../new_figs/modbb_ex_prop_pulse.pdf
