#!/usr/bin/gnuplot
# This file is covered by the LICENSE file in the root of this project.

set terminal pngcairo size 700,700 enhanced
set output 'distribution.png'
set title 'Uniform distribution on a sphere'

set view equal xyz
unset border
unset xtics
unset ytics
unset ztics

splot 'distribution.txt' with points pointtype 7 notitle
