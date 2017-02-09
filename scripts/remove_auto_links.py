import sys
import gzip

            

def remove_auto_linked_sequences(file_name):
    if "gz" in file_name:
        sequencefile=gzip.open(file_name,"r")
    else: 
        sequencefile=open(file_name,"r")
        
    #220:9174136-77-77.777779 4262888-97-97.979797 220-99-100.000000 
    
    for line in sequencefile.readlines():
        line = line.rstrip()
        if line[0]=='#':                    #Conserve the header file
            print (line)
            continue
        current_id=line.split(":")[0]       #220
        line=line.split(":")[1]             #9174136-77-77.777779 4262888-97-97.979797 220-99-100.000000
        tab_line = line.split()             #split the rest of the line
        # print (tab_line) #DEB
        if len(tab_line)==1: continue       # nothing to print, line is something like 220-99-100.000000
        toprint=current_id+":"
        for target in tab_line:
            if target.split("-")[0]==current_id: continue
            toprint+=target
        print (toprint)
    
    
if len(sys.argv)<2 :
    print ("USAGE")
    print (" Used to remove links between a sequence and itself in the SRC_linker format")
    print (" For other needs, contact pierre.peterlongo@inria.fr")
    print ("COMMAND")
    print (sys.argv[0],"<SRC_linker out file>")
    sys.exit(1)
remove_auto_linked_sequences(sys.argv[1])

