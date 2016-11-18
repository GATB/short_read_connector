import sys
import gzip
import canonical


def print_nodes_and_index_prefix_and_suffix(bank_file_name, len_prefix):

    prefixes={}
    suffixes={}
    if "gz" in bank_file_name:
        sequencefile=gzip.open(bank_file_name,"r")
    else:
        sequencefile=open(bank_file_name,"r")
        
    # i=1
    line=sequencefile.readline()
    if not line: 
        print ("Can't open file", bank_file_name)
        exit(1)
    if line[0]!='@' and line[0]!='>': 
        print ("File", bank_file_name, "not correctly formatted")
        exit(1)

    linesperread=2 #fasta by default
    if line[0]=='@': linesperread=4 # fastq
    
    sequencefile.seek(0)
    id_sequence=0
    while True:
        header=sequencefile.readline().rstrip() # read header
        if not header: break
        sequence = sequencefile.readline().rstrip() # read sequence
        prefixes[id_sequence]=sequence[0:len_prefix]
        suffixes[id_sequence]=sequence[-len_prefix-1:-1]
        # <node label="4" id="-2">
#           <att type="string" name="canonicalName" value="4"/>
#           <att type="string" name="comment" value="HELIUM_61LD0AAXX_61LD0AAXX:5:104:14649:7889#0"/>
#           <att type="integer" name="coverage" value="103"/>
#           <att type="integer" name="length" value="122"/>
#           <att type="string" name="sequence" value="TTCCTCTGATGAGTTGATGACGAAAATTAAAGATTCCACGCCGAATCTTTGTTTATTTTGCCAACCTAATAGTTATAAAAAATTCGCCGAGTTTCAGTCCTTTTAACAAAATCCTATTTAAC"/>
#           <att type="string" name="type" value="right"/>
#           <att type="string" name="vizmap:test NODE_FILL_COLOR" value="java.awt.Color[r=51,g=204,b=0]"/>
#           <graphics type="RECTANGLE" h="19.862258911132812" w="19.862258911132812" x="97.49642944335938" y="105.61238098144531" fill="#00cc00" width="3" outline="#000000" cy:nodeTransparency="1.0" cy:nodeLabelFont="Default-0-12" cy:borderLineType="solid"/>
#         </node>
        #<node id="0" label="AGAGGGATCTGTTATAGATCAAACTGCCATGACTCTTTTCCAGGGAGTATGGGTCTCCGCGTTAGTCTCGCTGAATAACTTTCGGTAAGACGCAAGCGGCCAGTATGCCAAA">
        #</node>
        print ("<node id=\""+str(id_sequence)+"\" label=\""+header+"\">")
        
        print ("<att type=\"string\" name=\"start_stop\" value=\""+sequence[0:10]+"..."+sequence[-11:-1]+"\"/>")
        print ("<att type=\"integer\" name=\"length\" value=\""+str(len(sequence))+"\"/>")
        print ("</node>")
        # print (str(id_sequence)+"[ label=\""+header+"_len_"+str(len(sequence))+"_"+sequence[0:10]+"..."+sequence[-11:-1]+"\"];")
        for i in range(linesperread-2): 
            line=sequencefile.readline() #don't care
        id_sequence+=1
    # print "max=",t
    sequencefile.close()
    return prefixes,suffixes



def convert_SRC_linker_output(remove_similar_reads, SRC_linker_output_file_name, prefixes,suffixes):
    #OPEN SRC OUTPUT
    if "gz" in SRC_linker_output_file_name:
        srcfile=gzip.open(SRC_linker_output_file_name,"r")
    else:
        srcfile=open(SRC_linker_output_file_name,"r")

    #787:787-100.000000 718-97.342522

    for line in srcfile.readlines():
        if line[0]=='#': #header
            continue
        line=line.rstrip()
        query_read_id=int(line.split(':')[0])-1 # -1 as the read ids are 1 based in SRC
        
        targets=line.split(':')[1].split(' ')
        for target in targets:
            target_read_id=int(target.split('-')[0])-1 # -1 as the read ids are 1 based in SRC
            # target_read_kmer_covers = int(target.split('-')[1])
            target_read_similarity=round(float(target.split('-')[1]),2)
            if target_read_id<query_read_id: continue # very specific to the assembly of perfect overlapping contigs
            if remove_similar_reads: 
                if target_read_id == query_read_id: continue

                if prefixes[target_read_id] == prefixes[query_read_id]: continue
                if suffixes[target_read_id ]== suffixes[query_read_id]: continue
                if canonical.rev_comp(prefixes[target_read_id]) == suffixes[query_read_id]: continue
                if canonical.rev_comp(suffixes[target_read_id]) == prefixes[query_read_id]: continue
            #print line
            #<edge source="0" target="1" label="ff">
            #</edge>
            prefix_suffix="False"
            # print (suffixes[target_read_id])
            # print (prefixes[query_read_id])
            if prefixes[target_read_id] == suffixes[query_read_id]: prefix_suffix="True"
            elif suffixes[target_read_id] == prefixes[query_read_id]: prefix_suffix="True"
            elif canonical.rev_comp(prefixes[target_read_id]) == prefixes[query_read_id]: prefix_suffix="True"
            elif canonical.rev_comp(suffixes[target_read_id]) == suffixes[query_read_id]: prefix_suffix="True"
            
            print ("<edge source=\""+str(query_read_id)+"\" target=\""+str(target_read_id)+"\" label=\""+str(target_read_similarity)+"\">")
            print ("<att type=\"boolean\" name=\"prefix_suffix\" value=\""+prefix_suffix+"\"/>")
            print ("</edge>")
            # print (str(query_read_id)+" -> "+str(target_read_id)+" [label=\""+str(target_read_similarity)+"\"];")
            # print "################################ NEW COUPLE ("+str(query_read_id)+"/"+str(target_read_id)+":"+str(target_read_similarity)+" percent similarity)+##############################"
    srcfile.close()

    
    
if len(sys.argv)<3 or len(sys.argv)>4:
     print ("USAGE")
     print (" Used to transform the SRC_linker output into a xgmml format")
     print (" Each couple of similar reads are written (starting with an horrible \"###############################\" line")
     print (" This tool must not be used when SRC_linker was used to compare two distinct read sets but only when SRC_linker was used to compare a read set against itself")
     print (" For other needs, contact pierre.peterlongo@inria.fr")
     print ("COMMAND")
     print ( sys.argv[0],"<banq file (fasta format, each sequence on ONE line, gzipped or not)> <SRC_linker output file> <optional 'R': if present similarity between read and itself are removed>")
     sys.exit(1)


read_bank=sys.argv[1]
SRC_output=sys.argv[2]
remove_similar_reads = False
if len(sys.argv) == 4:
    if sys.argv[3]=='R': remove_similar_reads = True
print ("<graph label=\"res_rconnector\"\nxmlns:dc=\"http://purl.org/dc/elements/1.1/\" \nxmlns:xlink=\"http://www.w3.org/1999/xlink\" \nxmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" \n xmlns:cy=\"http://www.cytoscape.org\" \n xmlns=\"http://www.cs.rpi.edu/XGMML\"  \ndirected=\"0\">")
prefixes,suffixes=print_nodes_and_index_prefix_and_suffix(read_bank,60)
convert_SRC_linker_output(remove_similar_reads, SRC_output, prefixes,suffixes)
print ("</graph>")