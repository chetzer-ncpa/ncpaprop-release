unset grid
set key top right opaque samplen 1 font ",10"
set style data lines

set ylabel "Altitude [km]" offset 1

set term pngcairo enh color font ",14" size 820,500 lw 2
set out "../ePape_ex2_profiles.png"

set multiplot
set size 0.5,1

set origin 0,0
set title "Easterly Wind Speed"
set xlabel "Wind speed [m/s]"
plot [*:*] [0:140] \
	dir."toy_profiles/toy_soundspeed_jetstr0.dat" using 2:1 title "0 km",\
	dir."toy_profiles/toy_soundspeed_jetstr2.dat" using 2:1 title "350 km",\
	dir."toy_profiles/toy_soundspeed_jetstr4.dat" using 2:1 title "500 km",\
	dir."toy_profiles/toy_soundspeed_jetstr6.dat" using 2:1 title "650 km",\
	dir."toy_profiles/toy_soundspeed_jetstr8.dat" using 2:1 title "800 km" lt 6

set origin 0.5,0
set title "Easterly Effective Sound Speed"
set xlabel "Sound speed [m/s]"
plot [250:450] [0:140] \
	dir."toy_profiles/toy_soundspeed_jetstr0.dat" using (sqrt(1.4*$7/$6/10) +$2):1 title "0 km",\
	dir."toy_profiles/toy_soundspeed_jetstr2.dat" using (sqrt(1.4*$7/$6/10) +$2):1 title "350 km",\
	dir."toy_profiles/toy_soundspeed_jetstr4.dat" using (sqrt(1.4*$7/$6/10) +$2):1 title "500 km",\
	dir."toy_profiles/toy_soundspeed_jetstr6.dat" using (sqrt(1.4*$7/$6/10) +$2):1 title "650 km",\
	dir."toy_profiles/toy_soundspeed_jetstr8.dat" using (sqrt(1.4*$7/$6/10) +$2):1 title "800 km" lt 6
unset multiplot

set origin 0,0
set size 1,1

