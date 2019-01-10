cmake_minimum_required(VERSION 3.10)
project(CForwardModel_VAM)

set(CMAKE_C_STANDARD 99)

set(SOURCE_FILES main.c forwardmodel.h initialization.c forwardmodel.c gnssr.h surface.c ddm.c antenna.c math.c wind.c geom.c specular.c coord.c grid.c debug.c saveFile.c cygnss.h cygnss.c)
add_executable(CForwardModel_VAM ${SOURCE_FILES})

target_link_libraries(CForwardModel_VAM fftw3)
target_link_libraries(CForwardModel_VAM netcdf)
