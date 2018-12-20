test1: simulation.o
	g++ simulation.o -o test1

simulation.o: simulation.cpp cdlasso.h function.h
	g++ -c simulation.cpp -std=c++11

real: real_data.o
	g++ real_data.o -o real

real_data.o: real_data.cpp cdlasso.h function.h
	g++ -c real_data.cpp -std=c++11

clean:
	rm *.o test1 real