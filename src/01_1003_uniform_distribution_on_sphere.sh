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
set output '01_1003_uniform_distribution_on_sphere.png'

set view equal xyz
unset border
unset xtics
unset ytics
unset ztics

splot '01_1003_uniform_distribution_on_sphere.txt' with points pointtype 7 notitle
