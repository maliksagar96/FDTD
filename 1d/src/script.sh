cd data
./clear_all.sh
cd ..
rm ./a.out
# g++ 01_fdtd_abc.cpp
# g++ 02_fdtd_naive_pml.cpp
# g++ 03_fdtd_graded_pml.cpp
g++ 04_fdtd_CMPL.cpp
# g++ 05_fdtd_CMPL_frequency.cpp
./a.out 
python3 field_plotter.py
vlc Ex_field_movie.mp4
