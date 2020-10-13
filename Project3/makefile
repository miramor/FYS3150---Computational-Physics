all: compile execute

compile:
	c++ -std=c++17 main.cpp solver.cpp planet.cpp -o main.exe -larmadillo -O3

execute:
	#./main.py systemA 1 E
	./main.exe systemA 1 E
	#./plot_results.py systemA 1 E
