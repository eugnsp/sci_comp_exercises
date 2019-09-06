#!/usr/bin/gnuplot

# *********************************************************************
# Uniform distribution on sphere
# ------------------------------
# Problems and solutions in scientific computing by W.-H. Steeb et al.
# Chapter 10, problem 3

# Generate uniform sampling points on a sphere.

# This file is covered by the LICENSE file in the root of this project.
# **********************************************************************

set terminal pngcairo size 700,700 enhanced
set output '01_1011_ising_model_lattice.png'

set xrange [1:100]
set yrange [1:100]
unset cbtics
unset colorbox
unset border
unset xtics
unset ytics

plot '01_1011_ising_model_lattice.txt' matrix with image
