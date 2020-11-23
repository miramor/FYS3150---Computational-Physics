all: compile

compile:
	g++ -std=c++11 -O3 main.cpp isingModel.cpp -o output.exe -fopenmp

execute:
	./output.exe

clean:
	rm *.o output.exe
