| **Linux** | **Mac OSX** |
|-----------|-------------|
|           |             |
[![Build Status](https://ci.inria.fr/gatb-core/view/RConnector/job/tool-rconnector-build-debian7-64bits-gcc-4.7/badge/icon)](https://ci.inria.fr/gatb-core/view/RConnector/job/tool-rconnector-build-debian7-64bits-gcc-4.7/) | [![Build Status](https://ci.inria.fr/gatb-core/view/RConnector/job/tool-rconnector-build-macos-10.9.5-gcc-4.2.1/badge/icon)](https://ci.inria.fr/gatb-core/view/RConnector/job/tool-rconnector-build-macos-10.9.5-gcc-4.2.1/)

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

What is Short Read Connector (SRC)?
===================================

Short read connector enables the comparisons of two read sets *B* and *Q*. For
each read from *Q* it provides either: 

* ***short_read_connector_counter***: The number of occurrences of each *k*-mers of the read in the set *B*,
   or
* ***short_read_connector_linker***: A list of reads from *B* that share enough *k*-mers with the (a window of) the tested read from *A*.


***short_read_connector_counter*** and ***short_read_connector_linker*** both output a text file, as described in the corresponding sections.

These text files can be downstream exploited as described section "Back to reads"

### Citation

Marchet, Camille, et al. "A resource-frugal probabilistic dictionary and applications in bioinformatics." *Discrete Applied Mathematics* 274 (2020): 92-102.

Getting the latest source code
==============================

Requirements
------------

* CMake 3.1.0+; see http://www.cmake.org/cmake/resources/software.html

* c++ compiler; compilation was tested with gcc and g++ version\>=4.5 (Linux) and
  clang version\>=4.1 (Mac OSX).

Instructions
------------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get a local copy of source code
git clone --recursive https://github.com/GATB/rconnector.git

# compile the code an run a simple test on your computer
cd rconnector
sh INSTALL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Getting a binary stable release
===============================

Binary release for Linux and Mac OSX are provided within the "Releases" tab on
Github/rconnector web page.

#short_read_connector_counter

short_read_connector_counter: quick start
-----------

Run a simple test counting for each reads from data/c2.fasta.gz, the number of occurrences of each of its kmers (k=31 by default) in data/c1.fasta.gz. Kmers indexed from data/c1.fasta.gz are those occurring at least 2 times (by default).

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# First index the kmers from data/c1.fasta.gz:
sh ../short_read_connector_counter.sh index -b data/c1.fasta.gz -i my_counter_index.dumped

# Once indexation is made once, multiple queries may be performed as: 
ls data/c2.fasta.gz > fof.txt # creates a file of files, check the section "Input read sets" for details
sh ../short_read_connector_counter.sh query -i my_counter_index.dumped -q fof.txt  

# Cat short_read_connector_res.txt to check the results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

short_read_connector_counter: options
-------

### Indexation options

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Usage: sh short_read_connector_counter.sh index -b read_file -i dumped_index_name [OPTIONS]
   -b read_files for bank
     Example: -b data/c1.fasta.gz
   -i <string>. File of the index file to be created. Example "my_index.dumped"
   -k <int> value (<32). Set the length of used kmers. Must fit the compiled value. Default=31
   -f <int> value. Fingerprint size. Size of the key associated to each indexed value, limiting false positives. Default=12
   -G <int> value. gamma value. MPHF expert users parameter - Default=2
   -a <int> kmer_abundance_min (kmer from bank seen less than this value both in the bank are not indexed). Default=2
   -l Keep low complexity regions (default false)
   -t <int> number of thread used. Default=0 (all)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Query options

```
Usage: sh short_read_connector_counter.sh query -i indexed_file -q read_file_of_files [OPTIONS]
MANDATORY:
 -i index created by short_read_connector_counter.sh index (read only)
   Example: -i my_index.dumped
 -q read_file_of_files for query
   Example: -q data/fof.txt (with fof being a file of file descriptor)
OPTIONS:
   -w <int> window_size. See option -s. If the windows size is zero (default value), then the full read is considered. Default=0
   -s <int> kmer_threshold: Minimal percentage of shared kmer span for considering a query read as similar to a data set.  
            The kmer span is the number of bases from the read query covered by a kmer shared with the bank.
   -l Keep low complexity regions (default false)
            If a read of length 80 has a kmer-span of 60, then the percentage of shared kmer span is 75%. If a least a windows (of size "windows_size") contains at least kmer_threshold percent of positions covered by shared kmers, the read is output.
   -p <string> prefix. All out files will start with this prefix. Default="short_read_connector_res"
   -t <int> number of thread used. Default=0 (all)
```

## short_read_connector_counter: output

Two first lines of the output file:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 #query_read_id (from bank data/c2.fasta.gz) mean median min max percentage_shared_positions -- number of shared 31mers with banq my_counter_index.dumped
0 3.614286 4 2 5 100.000000
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first line is the file header. The second line can be decomposed as:

-   0: id of the query read (from read set contained in fof.txt)

-   3.614286: mean number of occurrences of its k-mers (here with k=31) in the
    read set *data/c1.fasta.gz*

-   4: median number of occurrences of its k-mers (here with k=31) in the read
    set *data/c1.fasta.gz*

-   2: minimal number of occurrences of at least a kmer from read 0 in the read
    set *data/c1.fasta.gz*

-   5: maximal number of occurrences of at least a kmer from read 0 in the read
    set *data/c1fasta.gz*

-   100.0000: percentage of positions in the best window from the query read
    (here '0') covered by a kmer indexed in the bank.

Note: it is possible that `percentage_shared_positions` is equal to 100% while
`min` is equal to 0. This means for instance that at a position `i`, the kmer
starting at this position is not shared but `i` is covered by a kmer starting at
another position.









#short_read_connector_linker

short_read_connector_linker: quick start
-----------

Run a simple test counting for each reads from data/c2.fasta.gz, the number of occurrences of each of its kmers (k=31 by default) in data/c1.fasta.gz. Kmers indexed from data/c1.fasta.gz are those occurring at least 2 times (by default).

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# First index the kmers from data/c1.fasta.gz:
sh ../short_read_connector_linker.sh index -b data/c1.fasta.gz -i my_linker_index.dumped

# Once indexation is made once, multiple queries may be performed as: 
ls data/c2.fasta.gz > fof.txt # creates a file of files, check the section "Input read sets" for details
sh ../short_read_connector_linker.sh query -i my_linker_index.dumped -q fof.txt  

# Cat short_read_connector_res.txt to check the results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

short_read_connector_linker: options
-------

### Indexation options

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Usage: sh short_read_connector_linker.sh index -b read_file  -i dumped_index_name [OPTIONS]
   -b read_files for bank
     Example: -b data/c1.fasta.gz
   -i <string>. File of the index file to be created. Example "my_index.dumped"
   -k <int> value (<32). Set the length of used kmers. Must fit the compiled value. Default=31
   -f <int> value. Fingerprint size. Size of the key associated to each indexed value, limiting false positives. Default=12
   -G <int> value. gamma value. MPHF expert users parameter - Default=2
   -a <int> kmer_abundance_min (kmer from bank seen less than this value both in the bank are not indexed). Default=2
   -l Keep low complexity regions (default false)
   -t <int> number of thread used. Default=0 (all)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Query options

```
Usage: sh short_read_connector_linker.sh query -i indexed_file -q read_file_of_files [OPTIONS]
MANDATORY:
 -i index created by short_read_connector_linker.sh index (read only)
   Example: -i my_index.dumped
 -q read_file_of_files for query
   Example: -q data/fof.txt (with fof being a file of file descriptor)
OPTIONS:
   -w <int> window_size. See option -s. If the windows size is zero (default value), then the full read is considered. Default=0
   -s <int> kmer_threshold: Minimal percentage of shared kmer span for considering 2 reads as similar.  
            The kmer span is the number of bases from the read query covered by a kmer shared with the target read.
            If a read of length 80 has a kmer-span of 60 with another read (of unkonwn size) from the bank, then the percentage of shared kmer span is 75%. If a least a windows (of size "windows_size") contains at least kmer_threshold percent of position covered by shared kmers, the read couple is output.
                 TRICK: with kmer_threshold<=0 a single kmer is sufficient in the linker mode to link two reads. 
   -l Keep low complexity regions (default false)
   -p <string> prefix. All out files will start with this prefix. Default="short_read_connector_res"
   -t <int> number of thread used. Default=0 (all)
   -r do not output precision about pair of similar reads. Only ids of reads from queries similar to at least one read from bank are output.
```

## short_read_connector_linker: output

Four first lines of the output file:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#query_read_id [target_read_id-kmer_span (k=31)-kmer_span query_percentage]* or U (unvalid read, containing not only ACGT characters or low complexity read)
#Target read set: my_linker_index.dumped
#Query read set number data/c2.fasta.gz
0:675-93-93.000000 808-89-89.000000
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The three first lines are the file header. The last line can be decomposed as:

-   0: id of the query read

-   675-93-93.000000: a target read and its peaces of information:

    -   675: id of the targeted read

    -   93: kmer-span (number of position of read 0 that is covered by at least
        a solid kmer present in read 675)

    -   93.000000: kmer-span ratio wrt to read 0 length (here 100) or to the
        window size (see option -w)

-   808-89-89.000000: a second targeted read and its pieces of information (and
    so on).

Input read sets
===============

We use file of files format for query. The input read sets are
provided using a file of file(s). The file of file(s) contains on each line a
read file or another file of file(s). Let's look to a few usual cases (italic
strings indicate the composition of a file):

-   Case1: I've a unique read set composed of a unique read file
    (`reads.fq.gz`). `fof.txt` contains`reads.fq.gz`

-   Case2: I've a unique read set composed of a couple of read files
    (`reads_R1.fq.gz` and`reads_R2.fq.gz`). This may be the case in case of pair
    end sequencing. `fof.txt` contains`fof_reads.txt`:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 with fof_reads.txt:

 * reads_R1.fq.gz
 * reads_R2.fq.gz
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-   Case3: I've two read sets each composed of a unique read file:
    `reads1.fq.gz` and `reads2.fq.gz`:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fof.txt:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
reads1.fq.gz
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
reads2.fq.gz
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-   Case4: I've two read sets each composed two read files: `reads1_R1.fq.gz`
    and `reads1_R2.fq.gz` and `reads2_R1.fq.gz`and`reads2_R2.fq.gz`:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fof.txt:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fof_reads1.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fof_reads2.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

with`fof_reads1.txt`:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
reads1_R1.fq.gz
reads1_R2.fq.gz
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

with fof_reads2.txt: 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
reads1_R1.fq.gz
reads1_R2.fq.gz
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

and so on...

 

Back to reads 
======

Using the short read connector output, it is possible to select reads from the
query file that meets user defined requirements. It is done using two
executables: `generate_bv` and `extract_reads_from_bv`.

-   `generate_bv` extract query reads whose similarity with the bank is higher
    than a user defined threshold. Its takes as input a short read connector
    output (counter or linker), a percentage threshold and it generates a
    boolean vector (.bv file) containing a readable comment and a compressed 0/1
    vector. 0 designs reads the do not meet the similarity threshold and *vice
    versa.*

    -   With short read connector **linker** it contains all reads from a query
        bank that are similar (depending on -s and -w options) to **at least** a
        read from the bank

    -   With short read connector **counter** it contains all reads from a query
        bank that are similar (depending on -s and -w options) to the bank

-   `extract_reads_from_bv` takes as input a read set and a related .bv file and
    outputs the reads corresponding to ones in the .bv file.

 

Contact
=======

Contact: Pierre Peterlongo: pierre.peterlongo@inria.fr
