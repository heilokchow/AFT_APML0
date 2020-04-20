COMPILE_DEF	=

ifndef N_SIM
N_SIM = 100
endif
COMPILE_DEF += -DN_SIM=$(N_SIM)

ifndef P_SIM
P_SIM = 200
endif
COMPILE_DEF += -DP_SIM=$(P_SIM)

ifndef N_REP
N_REP = 100
endif
COMPILE_DEF += -DN_REP=$(N_REP)

ifndef N_THREAD
N_THREAD = $(shell grep -c ^processor /proc/cpuinfo)
endif
COMPILE_DEF += -DN_THREAD=$(N_THREAD)

ifndef CV_STRUCT
CV_STRUCT = 0
endif
COMPILE_DEF += -DCV_STRUCT=$(CV_STRUCT)

.PHONY: run_real run_simulation

simulation: simulation.o
	g++ simulation.o -o test -lpthread -O2

simulation.o: simulation.cpp cdlasso.h function.h
	g++ -c simulation.cpp $(COMPILE_DEF) -std=c++11 -lpthread -O2

real: real_data.o
	g++ real_data.o -o real -std=c++11 -lpthread -O2

real_data.o: real_data.cpp cdlasso.h function.h
	g++ -c real_data.cpp -std=c++11 -lpthread -O2

test: real_data.cpp cdlasso.h function.h
	g++ -o test real_data.cpp -DILLUSTRATION -std=c++11 -lpthread -O2

clean:
	-rm *.o test

run_real: real
	./real

run_simulation: simulation
	./simulation	
