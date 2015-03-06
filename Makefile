CC=g++

ace: ace.cpp 
	$(CC) -D MAXREADLEN=251 -c ace.cpp -o ace.o -fopenmp
	$(CC) ace.o -o ace -fopenmp -lpthread

