#!/usr/bin/gnuplot
# This file is covered by the LICENSE file in the root of this project.

set terminal pngcairo size 1000,700 enhanced
set output 'daubechies.png'

set title 'Daubechies lossy wavelet transform (J = 10, n = 2^J = 1024)'

set xlabel 'x'
set ylabel 'f'

plot 'daubechies.txt' using 2:3 with line linewidth 2 title 'Original function', \
     'daubechies.txt' using 2:4 with line linewidth 2 title '{/Symbol e} = .1', \
     'daubechies.txt' using 2:5 with line linewidth 2 title '{/Symbol e} = .5'

