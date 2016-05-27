import sys
            

def uniformizator(file_name):
    sequencefile=open(file_name,"r")
        
    #31:1246-3 479-3 1043-3 820-3 
    # will be changed into: 
    #31: 479-3 820-3 1043-3 1246-3 
    for line in sequencefile.readlines():
        if line[0]=='#': continue
        dic={}
        tab_line = line.rstrip().split(':')[1].split()
        for id_and_count in tab_line:
            dic[int(id_and_count.split('-')[0])]=int(id_and_count.split('-')[1])
        stringline=line.rstrip().split(':')[0]+':'
        for key in sorted(dic.keys()): 
            stringline+=str(key)+'-'+str(dic[key])+' '
        print stringline
        # print sorted(dic.get)
        # for key in sorted(dic, key=dic.get):
          # print key, dic[key]
    
    
if len(sys.argv)<2 :
    print "USAGE"
    print " Used to uniformize a SRC_linker output"
    print " For each line, the target read ids are put in a uniform order"
    print " The header is removed"
    print " For other needs, contact pierre.peterlongo@inria.fr"
    print "COMMAND"
    print  sys.argv[0],"SRC_output"
    sys.exit(1)
uniformizator(sys.argv[1])

