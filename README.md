ACE
===

ACE - Accurate Correction of Errors
-----------------------------------

ACE corrects substitution errors in single/paired-read Illumina archives
using a k-mer trie. It uses multiple cores when available, requiring OpenMP.
ACE has been tested on Linux, but should compile on any system in which g++
and OpenMP are available.

To compile the program:

  g++ -c ace.cpp -o ace.o -fopenmp
  g++ ace.o -o ace -fopenmp -lpthread

or simply enter "make" to use the Makefile supplied.

To run ACE:

  ace G InputFile(s)

where InputFile(s) can be either a FASTA or FASTQ file and G is the estimated genome 
length (in bp). ACE produces file(s) with .corrected extension.

(C) 2014, Siavash Sheikhizadeh
    Bioinformatics Group, Wageningen University
    e-mail: siavash.sheikhizadehanari@wur.nl

