cd data
./clear_all.sh
cd ..
g++ fdtd_2d.cpp
./a.out
python3 field_plotter.py