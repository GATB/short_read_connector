#Short Read Connector (SRC)



# What is SRC?

TODO
# Getting the latest source code

## Requirements

CMake 2.6+; see http://www.cmake.org/cmake/resources/software.html

c++ compiler; compilation was tested with gcc and g++ version>=4.5 (Linux) and clang version>=4.1 (Mac OSX).

## Instructions

    # get a local copy of DiscoSnp source code
    todo
    
    # compile the code an run a simple test on your computer
    cd commet_linked
    sh INSTALL

#Getting a binary stable release

todo

# Quick start

Run a simple test looking for reads from data/c2.fasta.gz that share at least 20 kmers (k=25) with data/c1.fasta.gz. Kmers indexed from data/c1.fasta.gz are those occurring at least 2 times. 

    sh short_read_connector.sh -b data/c1.fasta.gz -q data/c2.fasta.gz 


#User manual

Todo

#Contact

Contact: Pierre Peterlongo: pierre.peterlongo@inria.fr
