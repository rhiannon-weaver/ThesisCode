all: testNetworkGraph nat singleNetwork

singleNetwork: testing/runNATnetwork.cpp mcmcdrivers.o parameterupdaters.o helperfunctions.o graphlinks.o networkobjects.o datparobjects.o maintenance_functions.o distributionlist.o poissonEM.o mergesplitupdaters.o countupdaters.o
	g++ -g testing/runNATnetwork.cpp mcmcdrivers.o parameterupdaters.o graphlinks.o networkobjects.o datparobjects.o distributionlist.o helperfunctions.o maintenance_functions.o poissonEM.o mergesplitupdaters.o countupdaters.o -lgsl -lgslcblas -lm -o singleNetwork

nat: testing/nat.cpp mcmcdrivers.o parameterupdaters.o helperfunctions.o graphlinks.o networkobjects.o datparobjects.o maintenance_functions.o distributionlist.o poissonEM.o mergesplitupdaters.o countupdaters.o
	g++ -g testing/nat.cpp mcmcdrivers.o parameterupdaters.o graphlinks.o networkobjects.o datparobjects.o distributionlist.o helperfunctions.o maintenance_functions.o poissonEM.o mergesplitupdaters.o countupdaters.o -lgsl -lgslcblas -lm -o nat

testNetworkGraph: testing/testNetworkGraph.cpp mcmcdrivers.o parameterupdaters.o helperfunctions.o graphlinks.o networkobjects.o datparobjects.o maintenance_functions.o distributionlist.o poissonEM.o mergesplitupdaters.o countupdaters.o
	g++ -g testing/testNetworkGraph.cpp mcmcdrivers.o parameterupdaters.o graphlinks.o networkobjects.o datparobjects.o distributionlist.o helperfunctions.o maintenance_functions.o poissonEM.o mergesplitupdaters.o countupdaters.o -lgsl -lgslcblas -lm -o testNetworkGraph 

countupdaters.o: countupdaters.h mergesplitupdaters.h mcmcdrivers.h cppfiles/countupdaters.cpp
	g++ -g -c cppfiles/countupdaters.cpp -o countupdaters.o
mergesplitupdaters.o: mergesplitupdaters.h helperfunctions.h distributionlist.h cppfiles/mergesplitupdaters.cpp
	g++ -g -c cppfiles/mergesplitupdaters.cpp -o mergesplitupdaters.o 
graphlinks.o: graphlinks.h networkobjects.h distributionlist.h cppfiles/graphlinks.cpp
	g++ -g -c cppfiles/graphlinks.cpp -o graphlinks.o
mcmcdrivers.o: mcmcdrivers.h parameterupdaters.h cppfiles/mcmcdrivers.cpp
	g++ -g -c cppfiles/mcmcdrivers.cpp -o mcmcdrivers.o
networkobjects.o: networkobjects.h datparobjects.h cppfiles/networkobjects.cpp
	g++ -g -c cppfiles/networkobjects.cpp -o networkobjects.o
parameterupdaters.o: parameterupdaters.h helperfunctions.h distributionlist.h cppfiles/parameterupdaters.cpp 
	g++ -g -c cppfiles/parameterupdaters.cpp -o parameterupdaters.o 
datparobjects.o: datparobjects.h maintenance_functions.h cppfiles/datparobjects.cpp
	g++ -g -c cppfiles/datparobjects.cpp -o datparobjects.o
helperfunctions.o: helperfunctions.h datparobjects.h cppfiles/helperfunctions.cpp
	g++ -g -c cppfiles/helperfunctions.cpp -o helperfunctions.o
distributionlist.o: distributionlist.h cppfiles/distributionlist.cpp
	g++ -g -c cppfiles/distributionlist.cpp -o distributionlist.o
maintenance_functions.o: maintenance_functions.h cppfiles/maintenance_functions.cpp
	g++ -g -c cppfiles/maintenance_functions.cpp -o maintenance_functions.o
poissonEM.o: poissonEM.h cppfiles/poissonEM.cpp
	g++ -g -c cppfiles/poissonEM.cpp -o poissonEM.o
clean:
	rm -f *.o
	rm -f testNetworkGraph
	rm -f nat
	rm -f singleNetwork

