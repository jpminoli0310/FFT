# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5


CMakeFiles/fft.out.dir/library.mod.proxy: CMakeFiles/fft.out.dir/Library.f90.o.provides
CMakeFiles/fft.out.dir/Library.f90.o.provides.build:
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod library CMakeFiles/fft.out.dir/library.mod.stamp GNU
	$(CMAKE_COMMAND) -E touch CMakeFiles/fft.out.dir/Library.f90.o.provides.build
CMakeFiles/fft.out.dir/build: CMakeFiles/fft.out.dir/Library.f90.o.provides.build

CMakeFiles/fft.out.dir/prueba.f90.o.requires: CMakeFiles/fft.out.dir/library.mod.proxy
CMakeFiles/fft.out.dir/prueba.f90.o: CMakeFiles/fft.out.dir/library.mod.stamp
