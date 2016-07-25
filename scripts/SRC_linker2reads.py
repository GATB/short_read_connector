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


def convert_SRC_linker_output(read_offsets, remove_similar_reads, SRC_linker_output_file_name, bank_file_name):
    #OPEN SRC OUTPUT
    if "gz" in SRC_linker_output_file_name:
        srcfile=gzip.open(SRC_linker_output_file_name,"r")
    else: 
        srcfile=open(SRC_linker_output_file_name,"r")
    
    #OPEN BANK OUTPUT
    #We open two streams, one of them (bankfile) is the 'query' from SRC linker (31: in the following example). This stream has less random accesses
    if "gz" in bank_file_name:
        bankfile=gzip.open(bank_file_name,"r")
        bankfile2=gzip.open(bank_file_name,"r")
    else: 
        bankfile=open(bank_file_name,"r")
        bankfile2=open(bank_file_name,"r")
        
    
    
    #31:1246-30 479-63 1043-53 820-83 
    
    for line in srcfile.readlines():
        if line[0]=='#': #header
            continue
        line=line.rstrip()
        query_read_id=int(line.split(':')[0])-1 # -1 as the read ids are 1 based in SRC
        targets=line.split(':')[1].split(' ')
        for target in targets:
            target_read_id=int(target.split('-')[0])-1 # -1 as the read ids are 1 based in SRC
            target_read_similarity=int(target.split('-')[1])
            if remove_similar_reads and target_read_id == query_read_id: continue
            #print line
            print "################################ NEW COUPLE ("+str(query_read_id)+"/"+str(target_read_id)+":"+str(target_read_similarity)+" percent similarity)+##############################"
            print get_read(bankfile,read_offsets[query_read_id]),
            print get_read(bankfile2,read_offsets[target_read_id]),
    srcfile.close()
    bankfile.close()
    bankfile2.close()
    
    
if len(sys.argv)<3 or len(sys.argv)>4:
     print "USAGE"
     print " Used to transform the SRC_linker output into a format that retrieves the similar reads (instead of only read ids) "
     print " Each couple of similar reads are written (starting with an horrible \"###############################\" line"
     print " This tool must not be used when SRC_linker was used to compare two distinct read sets but only when SRC_linker was used to compare a read set against itself"
     print " For other needs, contact pierre.peterlongo@inria.fr"
     print "COMMAND"
     print  sys.argv[0],"<banq file (fasta format, each sequence on ONE line, gzipped or not)> <SRC_linker output file> <optional 'R': if present similarity between read and itself are removed>"
     sys.exit(1)


read_bank=sys.argv[1]
SRC_output=sys.argv[2]
remove_similar_reads = False
if len(sys.argv) == 4:
    if sys.argv[3]=='R': remove_similar_reads = True
read_offsets=index_bank_offsets(read_bank)
convert_SRC_linker_output(read_offsets, remove_similar_reads, SRC_output, read_bank)
