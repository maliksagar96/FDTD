# cd data
# ./clear_all.sh
# cd ..
cd build
rm fdtd_async
make -j6
clear 
make -j6
./fdtd_async
cd ..
python3 moviemaker.py 0
# python3 moviemaker.py 2
vlc Hz_animation.mp4