#!/usr/bin/gnuplot

set terminal pngcairo size 1000,700 
set output 'solution.png'
set title 'Boundary value problem solution'
set grid

set xlabel "x"
set ylabel "y"
splot 'solution.dat' binary matrix with lines title '-{/Symbol D}u = cos(2x) sin(2y)'
