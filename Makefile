# Sample size
n = 200

# Dimension
p = 10000

# Number of paramerters selected by stage 2 of APML0
maxit = 50

# Number of replications for simulation
rep = 1

# Initial seed
seed = 1

# LASSO tunning parameter for test used
lambda = 0.01

# Testing Run (Y) or Simulation (N)
np = Y

.PHONY: run
test1: simulation.o
	g++ simulation.o -o test1 -O2

simulation.o: simulation.cpp cdlasso.h function.h
	g++ -c simulation.cpp -std=c++11 -O2

real: real_data.o
	g++ real_data.o -o real -O2

real_data.o: real_data.cpp cdlasso.h function.h
	g++ -c real_data.cpp -std=c++11 -O2

clean:
	rm *.o test1 real

run: test1
	./test1 $(n) $(p) $(maxit) $(rep) $(seed) $(lambda) $(np)