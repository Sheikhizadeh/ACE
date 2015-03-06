ace: ace.cpp 
	g++ -D MAXREADLEN=??? -c ace.cpp -o ace.o -fopenmp
	g++ ace.o -o ace -fopenmp -lpthread
