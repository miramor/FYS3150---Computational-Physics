LIB = -larmadillo
CPPflags = g++ -std=c++17 -O3 #-fopenmp

output: main.o isingModel.o
	${CPPflags} main.o isingModel.o ${LIB} -o output.exe

main.o: main.cpp
		${CPPflags} -c main.cpp

isingModel.o: isingModel.cpp isingModel.hpp
	${CPPflags} -c isingModel.cpp

clean:
	rm *.o output.exe

execute:
	./output.exe

all: output execute
#compile:
#	g++ -fopenmp -o main.exe test.cpp
