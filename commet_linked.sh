#!/bin/bash

version="1.0.0"


EDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
platform='mac'
unamestr=`uname`
if [[ "$unamestr" == 'Linux' ]]; then
   platform='linux'
fi

dsk_bin=$EDIR/thirdparty/dsk/bin/macosx/dsk
if [[ $platform == 'linux' ]]; then
       dsk_bin=$EDIR/thirdparty/dsk/bin/linux/dsk
fi


function help {
echo "commet_linked.sh. Compare reads from two read sets (distinct or not)"
echo "Version "$version
echo "Usage: sh commet_linked.sh -b read_file_of_files -q read_file_of_files [OPTIONS]"
echo -e "\tMANDATORY:"
echo -e "\t\t -b read_file_of_files for bank"
echo -e "\t\t    Example: -b data/c1.fasta.gz"
echo -e "\t\t -q read_file_of_files for query"
echo -e "\t\t    Example: -q data/c2.fasta.gz"

echo -e "\tOPTIONS:"
echo -e "\t\t -p prefix. All out files will start with this prefix. Default=\"commet_linked_res\""
echo -e "\t\t -g: with this option, if a file of solid kmer exists with same prefix name and same k value, then it is re-used and not re-computed."
echo -e "\t\t -k value. Set the length of used kmers. Must fit the compiled value. Default=31"
echo -e "\t\t -f value. Fingerprint size. Size of the key associated to each indexed value, limiting false positives. Default=8"
echo -e "\t\t -a: kmer abundance min (kmer from bank seen less than this value are not indexed). Default=2"
echo -e "\t\t -s: minimal number of kmer shared by two reads to be considered as similar. Default=3"
echo -e "\t\t -t: number of thread used. Default=1"
echo -e "\t\t -d: use disk over RAM (slower)"
echo -e "\t\t -c: use commet_count"
}



bank_set=""
query_set=""
kmer_size=31
abundance_min=2
fingerprint_size=8
kmer_threshold=3
core_used=0
prefix="commet_linked_res"
remove=1
diskMode=0
countMode=0

#######################################################################
#################### GET OPTIONS                #######################
#######################################################################
while getopts "hgb:q:p:k:a:s:t:f:dc" opt; do
case $opt in

h)
help
exit
;;

d)
echo "use disk mode"
diskMode=1
;;

c)

echo "use commet count"
countMode=1
;;

f)
echo "use fingerprint size: $OPTARG" >&2
fingerprint_size=$OPTARG
;;

b)
echo "use bank read set: $OPTARG" >&2
bank_set=$OPTARG
;;

q)
echo "use query read set: $OPTARG" >&2
query_set=$OPTARG
;;

p)
echo "use prefix=$OPTARG" >&2
prefix=$OPTARG
;;

g)
echo "reuse solid precomputed solid kmers if exists"
remove=0
;;

k)
echo "use k=$OPTARG" >&2
kmer_size=$OPTARG
;;


a)
echo "use abundance_min=$OPTARG" >&2
abundance_min=$OPTARG
;;

s)
echo "use kmer_threshold=$OPTARG" >&2
kmer_threshold=$OPTARG
;;

t)
echo "use $OPTARG threads">&2
core_used=$OPTARG
;;

       #
       # u)
       # echo "use at most $OPTARG cores" >&2
       # option_cores_gatb="-nb-cores $OPTARG"
       # option_cores_post_analysis="-t $OPTARG"
       # ;;

\?)
echo "Invalid option: -$OPTARG" >&2
exit 1
;;

:)
echo "Option -$OPTARG requires an argument." >&2
exit 1
;;
esac
done
#######################################################################
#################### END GET OPTIONS            #######################
#######################################################################

if [ -z "${bank_set}" ]; then
	echo "You must provide a bank read set (-b)"
help
exit 1
fi

if [ -z "${query_set}" ]; then
	echo "You must provide a query read set (-q)"
help
exit 1
fi

out_dsk=${prefix}"_solid_kmers_k"${kmer_size}".h5"
result_file=${prefix}".txt"


if [ $remove -eq 1 ]; then
	rm -f ${out_dsk}
fi

# Count kmers using dsk if file absent
if [ ! -e ${out_dsk} ]; then
	echo "allo"
       cmd="${dsk_bin} -file ${bank_set},${query_set} -kmer-size ${kmer_size} -abundance-min ${abundance_min} -out ${out_dsk} -solidity-kind all"
       echo ${cmd}
       ${cmd}
fi

#unsorted_result_file=${result_file}"_unsorted"
# Compare read sets


# COMMET_LINKED_RAM
if [ $diskMode -eq 0 ]; then
	if [ $countMode -eq 0 ]; then
    	cmd="$EDIR/build/tools/commet_linked_ram/commet_linked_ram -graph ${out_dsk}  -bank ${bank_set} -query ${query_set} -out ${result_file} -kmer_threshold ${kmer_threshold} -fingerprint_size ${fingerprint_size} -core ${core_used}"
    else
		# COMMET_COUNT
       	cmd="$EDIR/build/tools/commet_count/commet_count -graph ${out_dsk}  -bank ${bank_set} -query ${query_set} -out ${result_file} -kmer_threshold ${kmer_threshold} -fingerprint_size ${fingerprint_size} -core ${core_used}"
       fi
else
	# COMMET_LINKED_DISK
	cmd="$EDIR/build/tools/commet_linked_disk/commet_linked_disk -graph ${out_dsk}  -bank ${bank_set} -query ${query_set} -out ${result_file} -kmer_threshold ${kmer_threshold} -fingerprint_size ${fingerprint_size} -core ${core_used}"
fi



echo ${cmd}
time ${cmd}

# sort results
#sort -n ${unsorted_result_file} > ${result_file}


#rm -f ${unsorted_result_file}


echo "***********************************"
echo "comment_linked finished"
echo "results in:"
echo "   "${result_file}
echo "Contact: pierre.peterlongo@inria.fr"
echo "***********************************"
