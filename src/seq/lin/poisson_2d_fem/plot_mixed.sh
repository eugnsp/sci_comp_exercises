#!/usr/bin/gnuplot

set terminal pngcairo size 1000,700 
set output 'mixed.png'
set title "Mixed FEM solution (P_1 + DP@_0^2 elements) of the BVP -{/Symbol D}u(x, y) = cos(2x) sin(2y)"

set grid
set xlabel "x"
set ylabel "y"

set xrange [0:3]
set yrange [0:3]
set zrange [-.15:.35]
set xyplane at -.15

splot 'mixed1.dat' with lines title 'u', \
	  'mixed2.dat' using 1:2:(-.15):($3/4):($4/4):(0) with vectors title '{/:Bold {/Symbol s}} = {/Symbol \321}u'


