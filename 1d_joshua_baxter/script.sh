cd data
./clear_all
cd ..
rm ./a.out
g++ fdtd_1d.cpp
./a.out 
python3 field_plotter.py
vlc Ex_field_movie.mp4&
sleep 10
pkill vlc
echo
