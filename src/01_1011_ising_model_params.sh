#!/usr/bin/gnuplot

set terminal pngcairo size 1000,700 enhanced
set output '01_1011_ising_model_params.png'

set multiplot layout 2,2 title 'The 2D Ising model (15 x 15 lattice, J = 1)'

set xlabel "T"
set ylabel "⟨E⟩"
set key at graph 0.4,0.95

plot '01_1011_ising_model_mt0.txt' using 1:2 with line linewidth 2 title 'H = 0', \
 	 '01_1011_ising_model_mt1.txt' using 1:2 with line linewidth 2 title 'H = 0.1', \
 	 '01_1011_ising_model_mt2.txt' using 1:2 with line linewidth 2 title 'H = 0.5'

set ylabel "⟨|M|⟩"
unset key

plot '01_1011_ising_model_mt0.txt' using 1:3 with line linewidth 2 title 'H = 0', \
 	 '01_1011_ising_model_mt1.txt' using 1:3 with line linewidth 2 title 'H = 0.1', \
 	 '01_1011_ising_model_mt2.txt' using 1:3 with line linewidth 2 title 'H = 0.5'

set ylabel "C_V"

plot '01_1011_ising_model_mt0.txt' using 1:4 with line linewidth 2 title 'H = 0', \
 	 '01_1011_ising_model_mt1.txt' using 1:4 with line linewidth 2 title 'H = 0.1', \
 	 '01_1011_ising_model_mt2.txt' using 1:4 with line linewidth 2 title 'H = 0.5'

set ylabel "{/Symbol c}"

plot '01_1011_ising_model_mt0.txt' using 1:5 with line linewidth 2 title 'H = 0', \
 	 '01_1011_ising_model_mt1.txt' using 1:5 with line linewidth 2 title 'H = 0.1', \
 	 '01_1011_ising_model_mt2.txt' using 1:5 with line linewidth 2 title 'H = 0.5'
	 