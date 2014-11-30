ACE
===

ACE - Accurate Correction of Errors
-----------------------------------

ACE corrects substitution errors in single/paired-read Illumina archives
using a k-mer trie. It uses multiple cores when available, requiring OpenMP.
ACE has should be be compiled on any Linux and Mac OS X system in which g++
and OpenMP are available.

To compile the program:

  g++ -c ace.cpp -o ace.o -fopenmp
  g++ ace.o -o ace -fopenmp -lpthread

or simply enter "make" to use the Makefile supplied.

To run ACE:

  ace G InputFile(s)

where G is the estimated genome length (in bp) and InputFile(s) can be either a FASTA or FASTQ file(s). ACE produces file(s) with .corrected before their original extension.

(C) 2014, Siavash Sheikhizadeh

    Bioinformatics Group, Wageningen University
    e-mail: siavash.sheikhizadehanari@wur.nl

