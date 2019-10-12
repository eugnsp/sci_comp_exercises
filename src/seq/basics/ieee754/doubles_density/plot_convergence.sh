#!/usr/bin/gnuplot
# This file is covered by the LICENSE file in the root of this project.

set terminal pngcairo size 1000,700
set output 'doubles_density.png'
set title 'Relative density of double floating-point numbers'
set grid

set logscale x 2
set xlabel "d"
set ylabel "(next(d) - d) / d"

plot 'doubles_density.txt' using 2:3 with lines linewidth 2 notitle
