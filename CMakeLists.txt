cmake_minimum_required(VERSION 3.0)
project(fft Fortran)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

add_executable(fft.out prueba.f90 Library.f90)
target_link_libraries(fft.out ${OpenMP_LIBS})

