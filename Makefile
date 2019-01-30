.PHONY: run_real run_simulation
simulation: simulation.o
	g++ simulation.o -o test1 -lpthread -O2

simulation.o: simulation.cpp cdlasso.h function.h
	g++ -c simulation.cpp -std=c++11 -lpthread -O2

real: real_data.o
	g++ real_data.o -o real -std=c++11 -lpthread -O2

real_data.o: real_data.cpp cdlasso.h function.h
	g++ -c real_data.cpp -std=c++11 -lpthread -O2

clean:
	rm *.o simulation real

run_real: real
	./real

run_simulation: simulation
	./simulation	
