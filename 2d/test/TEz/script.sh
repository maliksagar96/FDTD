#!/bin/bash

# Build the project
rm FDTD_2D
cd ../../build
make 
clear 
make
cd ../test/TEz
../../build/FDTD_2D control_file.scf
cd ../../data
python3 imagesc_plotter.py
vlc Ez_movie.mp4
# python3 sigma_plotter.py