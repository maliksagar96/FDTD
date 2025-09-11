cd ../build
rm read_step
make -j6
clear
make -j6
cd ../test
../build/read_step control_file.json