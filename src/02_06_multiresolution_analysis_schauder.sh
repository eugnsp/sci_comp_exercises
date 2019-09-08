#!/usr/bin/gnuplot

# *********************************************************************
# Multiresolution analysis
# ------------------------
# An introduction to scientific computing by I.Danaila et al.
# Chapter 6
#
# Implement 1D forward and backward wavelet transforms for Haar,
# Schauder and Daubechies wavelets.
#
# This file is covered by the LICENSE file in the root of this project.
# **********************************************************************

set terminal pngcairo size 1000,700 enhanced
set output '02_06_multiresolution_analysis_schauder.png'

set title 'Schauder lossy wavelet transform (J = 10, n = 2^J = 1024)'

set xlabel 'x'
set ylabel 'f'

plot '02_06_multiresolution_analysis_schauder.txt' using 1:2 with line linewidth 2 title 'original function', \
     '02_06_multiresolution_analysis_schauder.txt' using 1:3 with line linewidth 2 title '{/Symbol e} = .1', \
     '02_06_multiresolution_analysis_schauder.txt' using 1:4 with line linewidth 2 title '{/Symbol e} = .5'
    
