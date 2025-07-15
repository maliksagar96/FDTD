cd ../build/
make
clear
make
cd ../test/
../build/FDTD_1D control_file.scf
cd ../data
python3 field_plotter.py
vlc Ex_field_movie.mp4&
sleep 5
powershell.exe Stop-Process -Name vlc -Force
