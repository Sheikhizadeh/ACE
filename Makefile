ace: ace.cpp 
	g++ -c ace.cpp -o ace.o -fopenmp
	g++ ace.o -o ace -fopenmp -lpthread
