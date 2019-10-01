#!/bin/sh
# This file is covered by the LICENSE file in the root of this project.

gnuplot -persist <<-EOF
	set terminal pngcairo size 1280,720
	set title 'BVP solution using the SOR method'

	set grid
	set xlabel "x"
	set ylabel "y"
	set zrange [-.2:.4]

	n = `ls sor_*.dat | sort -V | tail -n 1 | sed 's/[^0-9]*//g'`
	do for [i=0:n] {
	    set output sprintf('sor_%d.png', i)
	    unset label 1
    	set label 1 sprintf('Iteration = %d', i + 1) at screen .1,.8
	    splot sprintf('sor_%d.dat', i) binary matrix with lines title '-{/Symbol D}u(x, y) = cos(2x) sin(2y)'
	}
EOF

ffmpeg -framerate 5 -i "sor_%d.png" -c:v libx264 -r 24 -pix_fmt yuv420p "sor_solution.mp4"
