rm -f CMakeCache.txt cmake_install.cmake 
rm -rf CMakeFiles/

export CC=/opt/gcc/8.3.0/bin/gcc
cmake .
make
