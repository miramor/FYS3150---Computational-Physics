all: compile

compile:
	g++ -std=c++17 main.cpp isingModel.cpp -o output.exe -fopenmp

execute:
	./output.exe

clean:
	rm *.o output.exe


#LIB = -larmadillo
#CPPflags = g++ -std=c++17 -O3 -fopenmp

#output: main.o isingModel.o
#	${CPPflags} main.o isingModel.o ${LIB} -o output.exe

#	main.o: main.cpp
#	${CPPflags} -c main.cpp

#	isingModel.o: isingModel.cpp isingModel.hpp
#	${CPPflags} -c isingModel.cpp

#	clean:
#	rm *.o output.exe
