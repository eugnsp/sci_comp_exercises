#!/usr/bin/gnuplot

set terminal pngcairo size 1000,700 
set output 'mixed_sy.png'
set title "Mixed FEM solution (P_1 + DP@_0^2 elements) of the BVP -{/Symbol D}u(x, y) = cos(2x) sin(2y)"

set grid
set xlabel "x"
set ylabel "y"

set xrange [0:3]
set yrange [0:3]

splot 'mixed_s.dat' using 1:2:4 with lines title '{/Symbol s}_y = {/Symbol \266}u/{/Symbol \266}y',\
