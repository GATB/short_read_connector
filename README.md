# What is Short Read Connector (SRC)?
Short read connector enables the comparisons of two read sets *B* and *Q*. For each read from *Q* it provides either:
 * The number of occurrences of each *k*-mers of the read in the set *B* (SRC_counter)
 or
 * A list of reads from *B* that share enough *k*-mers with the tested read from *B* (SRC_linker)


# Getting the latest source code

## Requirements

CMake 2.6+; see http://www.cmake.org/cmake/resources/software.html

c++ compiler; compilation was tested with gcc and g++ version>=4.5 (Linux) and clang version>=4.1 (Mac OSX).

## Instructions


    # get a local copy of DiscoSnp source code
    git clone --recursive https://github.com/GATB/rconnector.git
    
    # compile the code an run a simple test on your computer
    cd gatb-rconnector
    sh INSTALL

# Getting a binary stable release

todo

# Quick start

Run a simple test looking for reads from data/c2.fasta.gz that share at least 20 kmers (k=25) with data/c1.fasta.gz. Kmers indexed from data/c1.fasta.gz are those occurring at least 2 times. 

    sh short_read_connector.sh -b data/c1.fasta.gz -q data/c2.fasta.gz 


# Usage
## Mimimal call
Calling SRC_linker between read sets *bank* and *query*:
    sh short_read_connector.sh -b bank -q query

## Options
	 -p prefix. All out files will start with this prefix. Default="short_read_connector_res"
	 -g: with this option, if a file of solid kmer exists with same prefix name and same k value, then it is re-used and not re-computed.
	 -k value. Set the length of used kmers. Must fit the compiled value. Default=31
	 -f value. Fingerprint size. Size of the key associated to each indexed value, limiting false positives. Default=12
	 -G value. gamma value. MPHF expert users parameter - Default=2
	 -a: kmer abundance min (kmer from bank seen less than this value are not indexed). Default=2
	 -s: minimal number of kmer shared by two reads to be considered as similar. Default=3
	 -t: number of thread used. Default=0
	 -d: use disk over RAM (slower and no impact with -c option)
	 -c: use short_read_connector_counter (SRC_counter)
	 
	 
## Input read sets
We use file of files format. The input read sets are provided using a file of file(s). The file of file(s) contains on each line a read file or another file of file(s).
Let's look to a few usual cases (italic strings indicate the composition of a file):
* Case1: I've a unique read set composed of a unique read file (reads.fq.gz). 
   fof.txt:
   reads.fq.gz
* Case2: I've a unique read set composed of a couple of read files (reads_R1.fq.gz and reads_R2.fq.gz). This may be the case in case of pair end sequencing.
   fof.txt:
    fof_reads.txt:
   with fof_reads.txt:
    reads_R1.fq.gz
    reads_R2.fq.gz
* Case3: I've two read sets each composed of a unique read file: reads1.fq.gz and reads2.fq.gz:
   fof.txt:
    reads1.fq.gz
    reads2.fq.gz
* Case4:  I've two read sets each composed two read files: reads1_R1.fq.gz and reads1_R2.fq.gz and  reads2_R1.fq.gz and reads2_R2.fq.gz:
   fof.txt:
    fof_reads1.txt
    fof_reads2.txt
   with fof_reads1.txt:
    reads1_R1.fq.gz
    reads1_R2.fq.gz
   with fof_reads2.txt:
    reads2_R1.fq.gz
    reads2_R2.fq.gz
* and so on...

   
#Contact

Contact: Pierre Peterlongo: pierre.peterlongo@inria.fr
