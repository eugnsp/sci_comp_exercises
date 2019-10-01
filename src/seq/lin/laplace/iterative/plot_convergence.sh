#!/usr/bin/gnuplot

set terminal pngcairo size 1000,700 
set output 'convergence.png'
set title 'Convergence of iterative methods'
set grid

set xlabel "Iteration"
set ylabel "log_{10} ||Lu_i - b||_{/Symbol \245}"

plot 'convergence.txt' using 1:2 with lines linewidth 2 title 'Jacobi', \
	 'convergence.txt' using 1:3 with lines linewidth 2 title 'Gauss-Seidel', \
	 'convergence.txt' using 1:4 with lines linewidth 2 title 'Successive over-relaxation (SOR)', \
	 'convergence.txt' using 1:5 with lines linewidth 2 title 'Even-odd SOR with Chebyshev acceleration'
