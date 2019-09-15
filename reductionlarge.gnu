set tics font ",30"
set terminal png size 1920, 1080
set pm3d map
set pm3d scansforward
set palette defined (0 "blue", 1 "cyan", 2 "green", 3 "yellow", 4 "red")
set nokey

set bmargin at screen 0.07
set tmargin at screen 0.95
set lmargin at screen 0.07
set rmargin at screen 0.85

set logscale x
set format x "%.2e"
set xrange[3.0e+4:3.0e+6]
set xtics(3.0e+4,1.0e+5,3.0e+5,1.0e+6,3.0e+6)
set yrange[0:1]
set logscale zcb
set zrange[1e-6:1]

set output 'reduction/FULLl2fluxerror.png'
set multiplot
set format cb "%.2e"
splot 'reduction/interpprint.dat' u 1:2:3 with points pt 6 ps 7.0
splot 'reduction/fullresults.dat' u 1:2:3 with points palette pt 7 ps 6.0
unset multiplot

set output 'reduction/FULLl2valueserror.png'
set multiplot
set format cb "%.2e"
splot 'reduction/interpprint.dat' u 1:2:3 with points pt 6 ps 7.0
splot 'reduction/fullresults.dat' u 1:2:4 with points palette pt 7 ps 6.0
unset multiplot

set output 'reduction/FULLcfluxerror.png'
set multiplot
set format cb "%.2e"
splot 'reduction/interpprint.dat' u 1:2:3 with points pt 6 ps 7.0
splot 'reduction/fullresults.dat' u 1:2:5 with points palette pt 7 ps 6.0
unset multiplot

set output 'reduction/FULLl2cvalueserror.png'
set multiplot
set format cb "%.2e"
splot 'reduction/interpprint.dat' u 1:2:3 with points pt 6 ps 7.0
splot 'reduction/fullresults.dat' u 1:2:6 with points palette pt 7 ps 6.0
unset multiplot

unset logscale zcb
set format cb "%.3e"
set output 'reduction/FULLfulltime.png'
set multiplot
set format cb "%4.2f"
set zrange[0:2]
splot 'reduction/interpprint.dat' u 1:2:3 with points pt 6 ps 7.0
unset zrange
splot 'reduction/fullresults.dat' u 1:2:7 with points palette pt 7 ps 6.0
unset multiplot

set output 'reduction/FULLredtime.png'
set multiplot
set format cb "%4.5f"
set zrange[0:2]
splot 'reduction/interpprint.dat' u 1:2:3 with points pt 6 ps 7.0
unset zrange
splot 'reduction/fullresults.dat' u 1:2:8 with points palette pt 7 ps 6.0
unset multiplot

