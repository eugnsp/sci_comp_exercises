#!/usr/bin/gnuplot
# This file is covered by the LICENSE file in the root of this project.

set terminal pngcairo size 1000,700 enhanced
set output 'params.png'
set multiplot layout 2,2 title 'The 2D Ising model (15 x 15 lattice, J = 1)'

set xlabel "T"
set ylabel "⟨E⟩"
set key at graph 0.4,0.95

plot 'mt0.txt' using 1:2 with line linewidth 2 title 'H = 0', \
 	 'mt1.txt' using 1:2 with line linewidth 2 title 'H = 0.1', \
 	 'mt2.txt' using 1:2 with line linewidth 2 title 'H = 0.5'

set ylabel "⟨|M|⟩"
unset key

plot 'mt0.txt' using 1:3 with line linewidth 2 title 'H = 0', \
 	 'mt1.txt' using 1:3 with line linewidth 2 title 'H = 0.1', \
 	 'mt2.txt' using 1:3 with line linewidth 2 title 'H = 0.5'

set ylabel "C_V"

plot 'mt0.txt' using 1:4 with line linewidth 2 title 'H = 0', \
 	 'mt1.txt' using 1:4 with line linewidth 2 title 'H = 0.1', \
 	 'mt2.txt' using 1:4 with line linewidth 2 title 'H = 0.5'

set ylabel "{/Symbol c}"

plot 'mt0.txt' using 1:5 with line linewidth 2 title 'H = 0', \
 	 'mt1.txt' using 1:5 with line linewidth 2 title 'H = 0.1', \
 	 'mt2.txt' using 1:5 with line linewidth 2 title 'H = 0.5'
