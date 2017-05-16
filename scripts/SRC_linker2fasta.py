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


def index_bank_offsets(bank_file_name):
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


def convert_SRC_linker_output(read_offsets, SRC_linker_output_file_name, read_bank_queries_file_name):
    #OPEN SRC OUTPUT
    if "gz" in SRC_linker_output_file_name:
        srcfile=gzip.open(SRC_linker_output_file_name,"r")
    else: 
        srcfile=open(SRC_linker_output_file_name,"r")
    
    #OPEN BANK OUTPUT
    #We open two streams, one of them (bankfile) is the 'query' from SRC linker (31: in the following example). This stream has less random accesses
    if "gz" in read_bank_queries_file_name:
        bankfile=gzip.open(read_bank_queries_file_name,"r")
    else: 
        bankfile=open(read_bank_queries_file_name,"r")
        
    
    
    #787:787-4064-100.000000 718-3956-97.342522 or
    #787 (with the commet like option)
    
    for line in srcfile.readlines():
        if line[0]=='#': #header
            continue
        line=line.rstrip()
        query_read_id=int(line.split(':')[0])
        #~ query_read_id=int(line.split(':')[0])-1 # -1 as the read ids are 1 based in SRC // NO LONGER
        print get_read(bankfile,read_offsets[query_read_id]),
    srcfile.close()
    bankfile.close()
    
    
if len(sys.argv)<3 or len(sys.argv)>4:
     print "USAGE"
     print " Used to transform the SRC_linker output into a set of reads from the query (only reads for the query that share enough similarity with at least one read from the bank"
     print " For other needs, contact pierre.peterlongo@inria.fr"
     print "COMMAND"
     print  sys.argv[0],"<query set file (fasta format, each sequence on ONE line, gzipped or not)> <SRC_linker output file>"
     sys.exit(1)


read_bank_queries_file_name=sys.argv[1]
SRC_linker_output_file_name=sys.argv[2]
read_offsets=index_bank_offsets(read_bank_queries_file_name)
convert_SRC_linker_output(read_offsets, SRC_linker_output_file_name, read_bank_queries_file_name)
