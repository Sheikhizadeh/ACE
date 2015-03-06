CC=g++-4.9

ace: ace.cpp 
	$(CC) -D MAXREADLEN=100 -c ace.cpp -o ace.o -fopenmp
	$(CC) ace.o -o ace -fopenmp -lpthread

