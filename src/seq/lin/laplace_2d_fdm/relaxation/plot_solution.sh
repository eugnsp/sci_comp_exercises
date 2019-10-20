#!/usr/bin/gnuplot
# This file is covered by the LICENSE file in the root of this project.

set terminal pngcairo size 1000,700
set output 'solution.png'
set title "FDM solution of the BVP -{/Symbol D}u(x, y) = cos(2x) sin(2y)"

set grid
set xyplane at -.1
set xlabel "x"
set ylabel "y"

splot 'solution.dat' binary matrix with lines title 'u'

