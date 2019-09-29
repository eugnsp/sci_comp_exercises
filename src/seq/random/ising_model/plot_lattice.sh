#!/usr/bin/gnuplot
# This file is covered by the LICENSE file in the root of this project.

set terminal pngcairo size 700,700 enhanced
set output 'lattice.png'

set xrange [1:100]
set yrange [1:100]
unset cbtics
unset colorbox
unset border
unset xtics
unset ytics

plot 'lattice.txt' matrix with image
