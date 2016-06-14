import sys
import gzip

            
    

def filter(file_name, threshold):
    if "gz" in file_name:
        sequencefile=gzip.open(file_name,"r")
    else: 
        sequencefile=open(file_name,"r")
        
    #qseqid sseqid evalue pident length qlen slen
    #OM-RGC.v1.000000024_Eukaryota	OM-RGC.v1.008747545_Eukaryota	2.22e-04	96.970	33	30690	825
    previous_query="rien"
    previous_target="rien"
    for line in sequencefile.readlines():
        tab_line = line.rstrip().split()
        if tab_line[0]==previous_query and tab_line[1]==previous_target: continue
        previous_query=tab_line[0]
        previous_target=tab_line[1]
        # the total size of aligned sequence aligned=(float(tab_line[3])*float(tab_line[4]))
        aligned=(float(tab_line[3])*float(tab_line[4]))
        # divided by
        ratio = aligned/float(tab_line[5])
        if  ratio <threshold: continue
        print tab_line[0]+'\t'+tab_line[1]+'\t'+tab_line[2]+'\t'+str(ratio)
    
    
if len(sys.argv)<3 :
    print "USAGE"
    print " Used to remove blast alignments whose aligned nucleotides portion are bellow a threshold wrt size of sequences"
    print " For other needs, contact pierre.peterlongo@inria.fr"
    print "COMMAND"
    print  sys.argv[0],"<blast cvs file with -outfmt \"6 qseqid sseqid evalue pident length\"> <threshold>"
    sys.exit(1)
filter(sys.argv[1], float(sys.argv[2]))

