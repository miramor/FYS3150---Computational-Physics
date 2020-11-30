LIB = -larmadillo
CPPflags = g++ -std=c++11

output: main.o SIRS.o
	${CPPflags} main.o SIRS.o ${LIB} -o output.exe

main.o: main.cpp
	${CPPflags} -c main.cpp

SIRS.o: SIRS.cpp SIRS.hpp
	${CPPflags} -c SIRS.cpp

clean:
	rm *.o output.exe
