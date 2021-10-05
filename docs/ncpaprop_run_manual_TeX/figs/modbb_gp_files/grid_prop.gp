c_max=0.320
T_0(x)=x/c_max
R_0=220.0
dR=20.0

scaler=1.0e6/1.68554

set key samplen 2 bottom right
set size 1,0.8
set term post enh eps color solid 22
set out "example_prop_pulse_grid.eps"

set title "Propagated Waveforms 220 to 300 km"
set xlabel "Travel Time [s]"
set ylabel "Normalized Amplitude 10^{-6} [ ]"

plot "code_outputs/mywavf_grid.dat" every ::2000:0:3500:0 using ($2+T_0(R_0)):(scaler*$3) lw 3 title "220 km",\
     "code_outputs/mywavf_grid.dat" every ::1800:1:3400:1 using ($2+T_0(R_0+dR)):(scaler*$3) lw 3 title "240 km",\
     "code_outputs/mywavf_grid.dat" every ::1600:2:3300:2 using ($2+T_0(R_0+2*dR)):(scaler*$3) lw 3 title "260 km",\
     "code_outputs/mywavf_grid.dat" every ::1400:3:3100:3 using ($2+T_0(R_0+3*dR)):(scaler*$3) lw 3 title "280 km",\
     "code_outputs/mywavf_grid.dat" every ::1100:4:3100:4 using ($2+T_0(R_0+4*dR)):(scaler*$3) lw 3 lt 7 title "300 km"

! epstopdf example_prop_pulse_grid.eps;rm example_prop_pulse_grid.eps
! mv example_prop_pulse_grid.pdf ../modbb_ex_prop_pulse_grid.pdf
! cp ../modbb_ex_prop_pulse_grid.pdf ../new_figs/modbb_ex_prop_pulse_grid.pdf