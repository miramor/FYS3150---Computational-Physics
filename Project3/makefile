LIB = -larmadillo
CPPflags = g++ -std=c++17 -O3

output: main.o planet.o solver.o
	${CPPflags} main.o planet.o solver.o ${LIB} -o output.exe

main.o: main.cpp
	${CPPflags} -c main.cpp

planet.o: planet.cpp planet.hpp
	${CPPflags} -c planet.cpp

solver.o: solver.cpp solver.hpp
	${CPPflags} -c solver.cpp

clean:
	rm *.o output.exe *.txt
