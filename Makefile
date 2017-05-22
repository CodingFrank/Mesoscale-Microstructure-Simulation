all: main.cpp grid.hpp grid.cpp tesselation.cpp tesselation.hpp output.cpp output.hpp
	mpic++ -I. -O3 main.cpp grid.cpp tesselation.cpp output.cpp -o mc.out