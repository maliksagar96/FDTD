cd ../build
rm fdtd_2d
make -j6
clear
make -j6
cd ../test
../build/fdtd_2d control_file.json
python3 moviemaker.py 1
vlc Ey_animation.mp4