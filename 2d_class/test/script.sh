cd data/
./clear_all.sh
cd ..
cd ../build
rm fdtd_2d
make -j6
clear
make -j6
cd ../test
../build/fdtd_2d control_file.json
# python3 moviemaker.py 5
# vlc Hz_animation.mp4

python3 moviemaker.py 2
vlc Ez_animation.mp4

