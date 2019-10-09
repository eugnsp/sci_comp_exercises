#!/usr/bin/gnuplot
# This file is covered by the LICENSE file in the root of this project.

set terminal pngcairo size 1280,720
set output 'solution.png'
set title 'BVP solution'

set grid
set xlabel "x"
set ylabel "y"

splot 'solution.dat' binary matrix with lines title '-{/Symbol D}u(x, y) = cos(2x) sin(2y)'
