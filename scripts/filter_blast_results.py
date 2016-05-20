import sys
import gzip

            

def filter(file_name, threshold):
    if "gz" in file_name:
        sequencefile=gzip.open(file_name,"r")
    else: 
        sequencefile=open(file_name,"r")
        
    #OM-RGC.v1.000000024_Eukaryota	OM-RGC.v1.000000024_Eukaryota	1.67e-25	87.963	108
    
    for line in sequencefile.readlines():
        tab_line = line.rstrip().split()
        if (float(tab_line[-2]) * int(tab_line[-1])/float(100)) <threshold: continue
        print threshold, (float(tab_line[-2]) * int(tab_line[-1])/float(100)), line,
    
    
if len(sys.argv)<3 :
    print "USAGE"
    print " Used to remove blast whose aligned nucleotides are bellow a threshold"
    print " For other needs, contact pierre.peterlongo@inria.fr"
    print "COMMAND"
    print  sys.argv[0],"<blast cvs file with -outfmt \"6 qseqid sseqid evalue pident length\"> <threshold>"
    sys.exit(1)
filter(sys.argv[1], int(sys.argv[2]))

