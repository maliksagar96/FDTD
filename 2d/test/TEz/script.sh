#!/bin/bash

rm FDTD_2D
cd ../../build
make 
clear 
make
cd ../test/TEz
../../build/FDTD_2D control_file.scf
cd ../../data
python3 imagesc_plotter.py
vlc Ez_movie.mp4&
sleep 13
# pkill vlc
powershell.exe taskkill /IM vlc.exe /F
# python3 sigma_plotter.py