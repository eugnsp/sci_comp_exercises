#!/usr/bin/gnuplot

set terminal pngcairo size 700,700 enhanced
set output '01_1003_uniform_distribution_on_sphere.png'

set view equal xyz
unset border
unset xtics
unset ytics
unset ztics

splot '01_1003_uniform_distribution_on_sphere.txt' with points pointtype 7 notitle