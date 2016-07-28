| **Linux** | **Mac OSX** |
|-----------|-------------|
[![Build Status](https://ci.inria.fr/gatb-core/view/RConnector/job/tool-rconnector-build-debian7-64bits-gcc-4.7/badge/icon)](https://ci.inria.fr/gatb-core/view/RConnector/job/tool-rconnector-build-debian7-64bits-gcc-4.7/) | [![Build Status](https://ci.inria.fr/gatb-core/view/RConnector/job/tool-rconnector-build-macos-10.9.5-gcc-4.2.1/badge/icon)](https://ci.inria.fr/gatb-core/view/RConnector/job/tool-rconnector-build-macos-10.9.5-gcc-4.2.1/)

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)


# What is Short Read Connector (SRC)?
Short read connector enables the comparisons of two read sets *B* and *Q*. For each read from *Q* it provides either:
 * The number of occurrences of each *k*-mers of the read in the set *B* (SRC_counter)
 or
 * A list of reads from *B* that share enough *k*-mers with the tested read from *B* (SRC_linker)
 
**Citation** Camille Marchet, Antoine Limasset, Lucie Bittner, Pierre Peterlongo. A resource-frugal probabilistic
dictionary and applications in (meta)genomics. 2016. <hal-01322440>


# Getting the latest source code

## Requirements

CMake 2.6+; see http://www.cmake.org/cmake/resources/software.html

c++ compiler; compilation was tested with gcc and g++ version>=4.5 (Linux) and clang version>=4.1 (Mac OSX).

## Instructions


    # get a local copy of source code
    git clone --recursive https://github.com/GATB/rconnector.git
    
    # compile the code an run a simple test on your computer
    cd gatb-rconnector
    sh INSTALL

# Getting a binary stable release


Binary release for Linux and Mac OSX are provided within the "Releases" tab on Github/rconnector web page.

# Quick start

Run a simple test looking for reads from data/c2.fasta.gz that share at least 20 kmers (k=25) with data/c1.fasta.gz. Kmers indexed from data/c1.fasta.gz are those occurring at least 2 times. 

	 sh short_read_connector.sh -b data/c1.fasta.gz -q data/fof.txt


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
	 
## Output Format
### Short reads counter
Command:

	 sh short_read_connector.sh -b data/c1.fasta.gz -q data/fof.txt -c

Two first lines of the output file: 

	 #query_read_id mean median min max number of shared 31mers with banq read set data/c2.fasta.gz
	 0 3.614286 4 2 5
The first line is the file header. 
The second line can be decomposed as:
   * 0: id of the query read
   * 3.614286: mean number of occurrences of its k-mers (here with k=31) in the read set *data/c2.fasta.gz*
   * 4: median number of occurrences of its k-mers (here with k=31) in the read set *data/c2.fasta.gz*
   * 2: minimal number of occurrences of at least a kmer from read 0 in the read set *data/c2.fasta.gz*
   * 5: maximal number of occurrences of at least a kmer from read 0 in the read set *data/c2.fasta.gz*

### Short reads linker
Command:

	 sh short_read_connector.sh -b data/c1.fasta.gz -q data/fof.txt

Two first lines of the output file: 

	 #query_read_id [target_read_id-kmer_span (k=31)-kmer_span query percentage]* or U (unvalid read, containing not only ACGT characters or low complexity read)
	 1:676-93-93.000000 809-89-89.000000
The first line is the file header. 
The second line can be decomposed as:
   * 1: id of the query read
   * 676-93-93.000000: a target read and its peaces of information:
   	* 676: id of the targeted read
   	* 93: kmer-span (number of position of read 1 that is covered by at least a solid kmer present in read 676)
   	* 93.000000: kmer-span ratio wrt to read 1 length (here 100)
   * 809-89-89.000000: a second targeted read and its pieces of information (and so on)
   
   
## Input read sets
We use file of files format. The input read sets are provided using a file of file(s). The file of file(s) contains on each line a read file or another file of file(s).
Let's look to a few usual cases (italic strings indicate the composition of a file):
* Case1: I've a unique read set composed of a unique read file (reads.fq.gz). 
   * fof.txt:
   * reads.fq.gz
* Case2: I've a unique read set composed of a couple of read files (reads_R1.fq.gz and reads_R2.fq.gz). This may be the case in case of pair end sequencing.
   * fof.txt:
     * fof_reads.txt:
   
     with fof_reads.txt:
    
     * reads_R1.fq.gz
     * reads_R2.fq.gz
* Case3: I've two read sets each composed of a unique read file: reads1.fq.gz and reads2.fq.gz:
   * fof.txt:
     * reads1.fq.gz
     * reads2.fq.gz
* Case4:  I've two read sets each composed two read files: reads1_R1.fq.gz and reads1_R2.fq.gz and  reads2_R1.fq.gz and reads2_R2.fq.gz:
   * fof.txt:
     * fof_reads1.txt
     * fof_reads2.txt
  
   with fof_reads1.txt:
  
     * reads1_R1.fq.gz
     * reads1_R2.fq.gz
   
   with fof_reads2.txt:
     * reads2_R1.fq.gz
     * reads2_R2.fq.gz
* and so on...

   
#Contact

Contact: Pierre Peterlongo: pierre.peterlongo@inria.fr
