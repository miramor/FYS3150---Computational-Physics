
LIB = -larmadillo -llapack -lblas
CPPflags = c++ -std=c++17

output: main.o planet.o
${CPPflags} main.o planet.o ${LIB} -o output.exe

main.o: main.cpp
${CPPflags} -c main.cpp

planet.o: planet.cpp planet.hpp
${CPPflags} -c planet.cpp

solve.o: solver.cpp solver.hpp
${CPPflags} -c solve.cpp

clean:
rm *.o output.exe *.txt
