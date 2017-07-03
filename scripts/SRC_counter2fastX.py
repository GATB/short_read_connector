import sys
import gzip

def get_read(sequencefile,offset):
    sequencefile.seek(offset)
    read=""
    line=sequencefile.readline()
    if not line: 
        print "cannot read read at offset", offset
        exit(1)
    read+=line#include header
    read+=sequencefile.readline()#include sequence
    if read[0]=='>': return read
    if read[0]!='@': 
        print "read offset", offset, "does not start with @ or >"
        exit(1)
    read+=sequencefile.readline()#include header2
    read+=sequencefile.readline()#include quality
    return read


def index_offsets(bank_file_name):
    read_offsets= []

    if "gz" in bank_file_name:
        sequencefile=gzip.open(bank_file_name,"r")
    else: 
        sequencefile=open(bank_file_name,"r")
        
    # i=1
    line=sequencefile.readline()
    if not line: 
        print "Can't open file", bank_file_name
        exit(1)
    if line[0]!='@' and line[0]!='>': 
        print "File", bank_file_name, "not correctly formatted"
        exit(1)

    linesperread=2 #fasta by default
    if line[0]=='@': linesperread=4 # fastq
    
    sequencefile.seek(0)
    # t=0
    while True:
        offset=sequencefile.tell()
        line=sequencefile.readline()
        if not line: break
        # t+=1
        read_offsets.append(offset)
        for i in range(linesperread-1): line=sequencefile.readline()
    # print "max=",t
    sequencefile.close()
    return read_offsets


def convert_SRC_counter_output(SRC_linker_output_file_name, threshold):
    readfile=None
    srcfile = open(SRC_linker_output_file_name,"r")
    for line in srcfile.readlines():
        if line[0]=='#': # a header eg:#query_read_id (from bank ../data/c1.fasta.gz) mean median min max percentage_shared_positions -- number of shared 31mers with banq ../data/c1.fasta.gz
            read_file_name = line.split('(')[1].split(')')[0].split(" ")[2] # dirty, sorry.
            read_offsets=index_offsets(read_file_name)
            if readfile: readfile.close()
            if "gz" in read_file_name:
                readfile=gzip.open(read_file_name,"r")
            else: 
                readfile=open(read_file_name,"r")
        else:               # a counted line: 1 5.728571 6 0 10 84.000000: read id1 kmers occur in avg 5.72 times in the bank, mean is 6 times, min is 0 and max is 10. 84% of the read is covered by kmers indexed in the bank. We apply the threshold on this last field. 
            line=line.rstrip().split(" ")
            if float(line[-1])>=threshold:
                 print get_read(readfile,read_offsets[int(line[0])]),
                
   
    readfile.close()
    srcfile.close()
    
    
if len(sys.argv)!=3:
     print "USAGE"
     print " Used to transform the SRC_counter output into a set of reads from the query (only reads for the query that share enough similarity with at least one read from the bank)"
     print " For other needs, contact pierre.peterlongo@inria.fr"
     print "COMMAND"
     print  sys.argv[0],"<SRC_linker output file> <threshold>"
     sys.exit(1)


# read_bank_queries_file_name=sys.argv[1]
SRC_linker_output_file_name=sys.argv[1]
convert_SRC_counter_output(SRC_linker_output_file_name,int(sys.argv[2]))
