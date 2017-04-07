import sys
import gzip

def get_read(sequencefile,offset,message):
    sequencefile.seek(offset)
    read=""
    line=sequencefile.readline()
    if not line: 
        print "cannot read read at offset", offset
        exit(1)
    read+=line[:-1]+"_"+message+"\n"#include header and message
    read+=sequencefile.readline()#include sequence
    if read[0]=='>': return read
    if read[0]!='@': 
        print "read offset", offset, "does not start with @ or >"
        exit(1)
    read+=sequencefile.readline()#include header2
    read+=sequencefile.readline()#include quality

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


def convert_SRC_counter_output(read_offsets, SRC_counter_output_file_name, read_bank_queries_file_name, min_avg_coverage_threshold, max_avg_coverage_threshold):
    #OPEN SRC OUTPUT
    if "gz" in SRC_counter_output_file_name:
        srcfile=gzip.open(SRC_counter_output_file_name,"r")
    else: 
        srcfile=open(SRC_counter_output_file_name,"r")
    
    #OPEN BANK OUTPUT
    #We open two streams, one of them (bankfile) is the 'query' from SRC counter (31: in the following example). This stream has less random accesses
    if "gz" in read_bank_queries_file_name:
        bankfile=gzip.open(read_bank_queries_file_name,"r")
    else: 
        bankfile=open(read_bank_queries_file_name,"r")
        
    
    
    #0 3.614286 4 2 5
    #id mean median min max

    
    for line in srcfile.readlines():
        if line[0]=='#': #header
            continue
        line=line.rstrip()
        
        avg_coverage=float(line.split()[1])
        if avg_coverage>=min_avg_coverage_threshold and avg_coverage<=max_avg_coverage_threshold:
            query_read_id=int(line.split()[0])
            print get_read(bankfile,read_offsets[query_read_id],"cov"+line[line.index(" "):].replace(' ', '_')),
    srcfile.close()
    bankfile.close()
    
    
if len(sys.argv)!=4 and len(sys.argv)!=5:
     print "USAGE"
     print " Used to transform the SRC_counter output into a set of reads from the query (only reads for the query that are covered enough."
     print " A minimal coverage threshold is provided and reads whose mean coverage is strictly lower than this threshold are not output."
     print " Optionally, user may provide a maximal coverage threshold. In this case, reads whose mean coverage is strictly higher than this threshold are not output."
     print " For other needs, contact pierre.peterlongo@inria.fr"
     print "COMMAND"
     print  sys.argv[0],"<query set file (fasta format, each sequence on ONE line, gzipped or not)> <SRC_counter output file> <min coverage threshold> [<max coverage threshold>]"
     sys.exit(1)


read_bank_queries_file_name=sys.argv[1]
SRC_counter_output_file_name=sys.argv[2]
min_avg_coverage_threshold=float(sys.argv[3])
max_avg_coverage_threshold=sys.maxint
if len(sys.argv)==5: max_avg_coverage_threshold=float(sys.argv[4])
read_offsets=index_bank_offsets(read_bank_queries_file_name)
convert_SRC_counter_output(read_offsets, SRC_counter_output_file_name, read_bank_queries_file_name,min_avg_coverage_threshold, max_avg_coverage_threshold)
