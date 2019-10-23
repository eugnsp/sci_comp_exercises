#!/usr/bin/gnuplot

set terminal pngcairo size 1000,700 
set output 'std.png'
set title "Standard FEM solution (P_4 elements) of the BVP -{/Symbol D}u(x, y) = cos(2x) sin(2y)"

set grid
set xlabel "x"
set ylabel "y"

set xrange [0:3]
set yrange [0:3]
set zrange [-.15:.35]
set xyplane at -.15

set multiplot
splot 'std_u.dat' with lines linecolor 1 notitle

set title ""
unset xlabel
unset ylabel
unset grid
unset tics
unset border

set dgrid3d 50,50,4
splot 'std_u_interp.dat' with lines linecolor 3 title 'u'

