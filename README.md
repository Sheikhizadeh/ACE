ACE
===

ACE - Accurate Correction of Errors
-----------------------------------

ACE corrects substitution errors in an Illumina archive using a k-mer trie. It uses multiple cores when available, requiring OpenMP. ACE has should be be compiled on any Linux and Mac OS X system in which g++ and OpenMP are available.It receives the genome length and the name of input and output files.
To compile the program run the provided makefile or type:

g++ -D MAXREADLEN=??? -c ace.cpp -o ace.o -fopenmp

Note that you should provide compiler with the maximum length of reads (MAXREADLEN=???).
Then type:

g++ ace.o -o ace -fopenmp -lpthread

To run ACE type:

./ace Genome_length(bp) Inputfile Outputfile

DOI: 10.1093/bioinformatics/btv332
(C) 2014, Siavash Sheikhizadeh

    Bioinformatics Group, Wageningen University
    e-mail: siavash.sheikhizadehanari@wur.nl
