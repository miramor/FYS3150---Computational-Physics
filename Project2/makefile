LIB = -larmadillo -llapack -lblas
CPPflags= c++ -
output: main.o jacobiAlgo.o
	${CPPflags} main.o jacobiAlgo.o ${LIB} -o output.exe

main.o: main.cpp
	${CPPflags} -c main.cpp

jacobiAlgo.o: jacobiAlgo.cpp jacobiAlgo.hpp
	${CPPflags} -c jacobiAlgo.cpp

clean:
	rm *.o output.exe *.txt
