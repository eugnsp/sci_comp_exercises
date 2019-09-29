#!/usr/bin/gnuplot
# This file is covered by the LICENSE file in the root of this project.

set terminal pngcairo size 1000,700 enhanced
set output 'schauder.png'

set title 'Schauder lossy wavelet transform (J = 10, n = 2^J = 1024)'

set xlabel 'x'
set ylabel 'f'

plot 'schauder.txt' using 2:3 with line linewidth 2 title 'Original function', \
     'schauder.txt' using 2:4 with line linewidth 2 title '{/Symbol e} = .1', \
     'schauder.txt' using 2:5 with line linewidth 2 title '{/Symbol e} = .5'

