import sys
import gzip

def index_headers(bank_file_name):
    headers= []
    headers.append("NULL") # reads are 1-indexed we add this dummy value
    if "gz" in bank_file_name:
        sequencefile=gzip.open(bank_file_name,"r")
    else: 
        sequencefile=open(bank_file_name,"r")
        
    # i=1
    for line in sequencefile.readlines():
        if line[0]=='>':
            # print line[1:],
            headers.append(line.rstrip()[1:])            #
            # print i
            # print line[1:],
            # print headers[i]
            # i+=1
    return headers
    


def sequence_sizes(bank_file_name):
    sizes= []
    sizes.append(0) # reads are 1-indexed we add this dummy value
    
    if "gz" in bank_file_name:
        sequencefile=gzip.open(bank_file_name,"r")
    else: 
        sequencefile=open(bank_file_name,"r")
        
    for line in sequencefile.readlines():
        if line[0]!='>':
            # print line[1:],
            sizes.append(len(line))
    return sizes
            

def convert_SRC_linker_output(headers, sizes, k, threshold, remove_similar_reads, SRC_linker_output_file_name):
    print "qseqid\tsseqid\tevalue\tpident"
    if "gz" in SRC_linker_output_file_name:
        sequencefile=gzip.open(SRC_linker_output_file_name,"r")
    else: 
        sequencefile=open(SRC_linker_output_file_name,"r")
        
    #31:1246-3 479-3 1043-3 820-3 
    
    for line in sequencefile.readlines():
        if line[0]=='#': #header
            continue
        line=line.rstrip()
        query_read_id=int(line.split(':')[0])
        targets=line.split(':')[1].split(' ')
        for target in targets:
            target_read_id=int(target.split('-')[0])
            if remove_similar_reads and target_read_id == query_read_id: continue
            coverage = 100*k*int(target.split('-')[1])/float(min(sizes[query_read_id],sizes[target_read_id]))
            print headers[query_read_id]+"\t"+headers[target_read_id]+"\t0.0\t%.4g"%(coverage)
    
    
if len(sys.argv)<5 or len(sys.argv)>6:
    print "USAGE"
    print " Used to transform the SRC_linker output into a format usable by mosaic finder"
    print " This tool must not be used when SRC_linker was used to compare two distinct read sets but only when SRC_linker was used to compare a read set against itself"
    print " For other needs, contact pierre.peterlongo@inria.fr"
    print "COMMAND"
    print  sys.argv[0],"<banq file (fasta format, each sequence on ONE line, gzipped or not)> <SRC_linker output file> <kmer size (used by SRC_linker)> <threshold: bellow this percentage, alignements whose coverage ratio is lower are not output> <optional 'R': if present similarity between read and itself are removed>"
    sys.exit(1)
headers=index_headers(sys.argv[1])
sizes=sequence_sizes(sys.argv[1])
k=int(sys.argv[3])
threshold=int(sys.argv[4])
remove_similar_reads = False
if len(sys.argv) == 6: 
    if sys.argv[5]=='R': 
        remove_similar_reads = True


convert_SRC_linker_output(headers,sizes,k, threshold, remove_similar_reads, sys.argv[2])

