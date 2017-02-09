import sys
import gzip

            

def remove_auto_linked_sequences(file_name):
    if "gz" in file_name:
        sequencefile=gzip.open(file_name,"r")
    else: 
        sequencefile=open(file_name,"r")
        
    #read18_contig0_position4295_M0_I0_D0_NG0______er0.01__indel0__rgeom0    read687_contig0_position4286_M1_I0_D0_NG0______er0.01__indel0__rgeom0   0.0     74.26
    
    for line in sequencefile.readlines():
        if line[0]=='#': continue
        tab_line = line.rstrip().split()
        if tab_line[0]==tab_line[1]: continue
        print line,
    
    
if len(sys.argv)<2 :
    print "USAGE"
    print " Used to remove links between a sequence and itself in the cvs format given to MosaicFinder"
    print " For other needs, contact pierre.peterlongo@inria.fr"
    print "COMMAND"
    print  sys.argv[0],"<cvs file>"
    sys.exit(1)
remove_auto_linked_sequences(sys.argv[1])

