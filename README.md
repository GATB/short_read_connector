| **Linux**                                                                                         | **Mac OSX**                                                                                       |
|---------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------|
|                                                                                                   |                                                                                                   |
| <https://ci.inria.fr/gatb-core/view/RConnector/job/tool-rconnector-build-debian7-64bits-gcc-4.7/> | <https://ci.inria.fr/gatb-core/view/RConnector/job/tool-rconnector-build-macos-10.9.5-gcc-4.2.1/> |

<http://www.gnu.org/licenses/agpl-3.0.en.html>

What is Short Read Connector (SRC)?
===================================

Short read connector enables the comparisons of two read sets *B* and *Q*. For
each read from *Q* it provides either: \* The number of occurrences of each
*k*-mers of the read in the set *B* (SRC_counter) or \* A list of reads from *B*
that share enough *k*-mers with the (a window of) the tested read from *A*
(SRC_linker)

**Citation** Camille Marchet, Antoine Limasset, Lucie Bittner, Pierre
Peterlongo. A resource-frugal probabilistic dictionary and applications in
(meta)genomics. 2016.

Getting the latest source code
==============================

Requirements
------------

CMake 2.6+; see http://www.cmake.org/cmake/resources/software.html

c++ compiler; compilation was tested with gcc and g++ version\>=4.5 (Linux) and
clang version\>=4.1 (Mac OSX).

Instructions
------------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get a local copy of source code
git clone --recursive https://github.com/GATB/rconnector.git

# compile the code an run a simple test on your computer
cd gatb-rconnector
sh INSTALL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Getting a binary stable release
===============================

Binary release for Linux and Mac OSX are provided within the "Releases" tab on
Github/rconnector web page.

Quick start
===========

Run a simple test looking for reads from data/c2.fasta.gz that share at least 20
kmers (k=25) with data/c1.fasta.gz. Kmers indexed from data/c1.fasta.gz are
those occurring at least 2 times.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 sh short_read_connector.sh -b data/c1.fasta.gz -q data/fof.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Usage
=====

Mimimal call
------------

Calling SRC_linker between read sets *bank* and *query*:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sh short_read_connector.sh -b bank -q query
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Options
-------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Usage: sh short_read_connector.sh -b read_file -q read_file_of_files [OPTIONS]
MANDATORY:
 -b read_files for bank
   Example: -b data/c1.fasta.gz
 -q read_file_of_files for query
   Example: -q data/fof.txt (with fof being a file of file descriptor)
OPTIONS:
   -c use short_read_connector_counter (SRC_counter)
 OPTIONS USING SRC_counter (-c) OR SRC_linker
   -w <int> window_size. See option -s. If the windows size is zero (default value), then the full read is considered. Default=0
   -p <string> prefix. All out files will start with this prefix. Default="short_read_connector_res"
   -g with this option, if a file of solid kmer exists with same prefix name and same k value, then it is re-used and not re-computed.
   -k <int> value. Set the length of used kmers. Must fit the compiled value. Default=31
   -f <int> value. Fingerprint size. Size of the key associated to each indexed value, limiting false positives. Default=12
   -G <int> value. gamma value. MPHF expert users parameter - Default=2
   -a <int> kmer_abundance_min (kmer from bank seen less than this value both in the bank are not indexed). Default=2
   -l Keep low complexity regions (default false)
   -t <int> number of thread used. Default=0 (all)
 OPTIONS USING SRC_counter (-c)
   -s <int> kmer_threshold: Minimal percentage of shared kmer span for considering a query read as similar to a data set.  
            The kmer span is the number of bases from the read query covered by a kmer shared with the bank.
            If a read of length 80 has a kmer-span of 60, then the percentage of shared kmer span is 75%. If a least a windows (of size "windows_size") contains at least kmer_threshold percent of positions covered by shared kmers, the read is output.
 OPTIONS USING SRC_linker (without -c)
   -s <int> kmer_threshold: Minimal percentage of shared kmer span for considering 2 reads as similar.  
            The kmer span is the number of bases from the read query covered by a kmer shared with the target read.
            If a read of length 80 has a kmer-span of 60 with another read (of unkonwn size) from the bank, then the percentage of shared kmer span is 75%. If a least a windows (of size "windows_size") contains at least kmer_threshold percent of position covered by shared kmers, the read couple is output.
   -A index kmers present at least 'kmer_abundance_min' times in the bank AND in the queries.
   -r (incompatible with SRC_counter), do not output precision about pair of similar reads. Only ids of reads from queries similar to at least one read from bank are output.
 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Input read sets
===============

Short read connector output is composed of two files:

-   A boolean vector file (.bv)

    -   With short read connector **linker** it contains all reads from a query
        bank that are similar (depending on -s and -w options) to at least a
        read from the bank

    -   With short read connector **counter** it contains all reads from a query
        bank that are similar (depending on -s and -w options) to the bank

    In order to obtain the corresponding reads, use the
    `build/bin/extract_reads_from_bv` executable

-   A text file, as described in the following section.

Output Format
-------------

### Short reads counter (with -c option)

Command:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 sh short_read_connector.sh -b data/c1.fasta.gz -q data/fof.txt -c
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two first lines of the output file:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 #query_read_id (from bank data/c2.fasta.gz) mean median min max percentage_shared_positions -- number of shared 31mers with banq data/c1.fasta.gz
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

 

### Short reads linker (without the -c option)

Command:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 sh short_read_connector.sh -b data/c1.fasta.gz -q data/fof.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Four first lines of the output file:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#query_read_id [target_read_id-kmer_span (k=31)-kmer_span query_percentage]* or U (unvalid read, containing not only ACGT characters or low complexity read)
#Target read set: data/c1.fasta.gz
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

Note that with the -r option, only the id of the queried and shared read is
output. In this example the line would be limited to

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 #query_read_id 
 0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 

Input read sets
---------------

We use file of files format for query and banks. The input read sets are
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

with fof_reads2.txt: \* reads2_R1.fq.gz \* reads2_R2.fq.gz \* and so on...

 

 

Contact
=======

Contact: Pierre Peterlongo: pierre.peterlongo\@inria.fr
