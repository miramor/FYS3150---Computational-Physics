all: compile execute

compile:
	g++ -std=c++11 -O3 main.cpp SIRS.cpp -larmadillo -o output.exe

execute:
	./output.exe

clean:
	rm *.o output.exe
