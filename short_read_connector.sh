#!/bin/bash

version="1.0.0"



EDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )
platform='mac'
unamestr=`uname`
if [[ "$unamestr" == 'Linux' ]]; then
    platform='linux'
fi

dsk_bin=$EDIR/thirdparty/dsk/bin/macosx/dsk
if [[ $platform == 'linux' ]]; then
    dsk_bin=$EDIR/thirdparty/dsk/bin/linux/dsk
fi
chmod a+x ${dsk_bin}


#BIN DIR:
if [ -d "$EDIR/build/" ] ; then # VERSION SOURCE COMPILED
    BIN_DIR=$EDIR/build/bin
else # VERSION BINARY
    BIN_DIR=$EDIR/bin
fi


function help {
    echo "short_read_connector.sh - Compare reads from two read sets (distinct or not)"
    echo "Version "$version
    echo "Usage: sh short_read_connector.sh -b read_file -q read_file_of_files [OPTIONS]"
    echo  "MANDATORY:"
    echo  " -b read_files for bank"
    echo  "   Example: -b data/c1.fasta.gz"
    echo  " -q read_file_of_files for query"
    echo  "   Example: -q data/fof.txt (with fof being a file of file descriptor)"

    echo  "OPTIONS:"
    echo  "   -c use short_read_connector_counter (SRC_counter)"
    echo  " OPTIONS USING SRC_counter (-c) OR SRC_linker"
    echo  "   -w <int> window_size. See option -s. If the windows size is zero (default value), then the full read is considered. Default=0"
    echo  "   -p <string> prefix. All out files will start with this prefix. Default=\"short_read_connector_res\""
    echo  "   -g with this option, if a file of solid kmer exists with same prefix name and same k value, then it is re-used and not re-computed."
    echo  "   -k <int> value. Set the length of used kmers. Must fit the compiled value. Default=31"
    echo  "   -f <int> value. Fingerprint size. Size of the key associated to each indexed value, limiting false positives. Default=12"
    echo  "   -G <int> value. gamma value. MPHF expert users parameter - Default=2"
    echo  "   -a <int> kmer_abundance_min (kmer from bank seen less than this value both in the bank are not indexed). Default=2"

    echo  "   -l Keep low complexity regions (default false)"
    echo  "   -t <int> number of thread used. Default=0 (all)"
#    echo  "   -d: use disk over RAM (slower and no impact with -c option)"

    echo  " OPTIONS USING SRC_counter (-c)"
    echo  "   -s <int> kmer_threshold: Minimal percentage of shared kmer span for considering a query read as similar to a data set.  "
    echo  "            The kmer span is the number of bases from the read query covered by a kmer shared with the bank."
    echo  "            If a read of length 80 has a kmer-span of 60, then the percentage of shared kmer span is 75%. If a least a windows (of size \"windows_size\") contains at least kmer_threshold percent of positions covered by shared kmers, the read is output."
    echo  " OPTIONS USING SRC_linker (without -c)"
    echo  "   -s <int> kmer_threshold: Minimal percentage of shared kmer span for considering 2 reads as similar.  "
    echo  "            The kmer span is the number of bases from the read query covered by a kmer shared with the target read."
    echo  "            If a read of length 80 has a kmer-span of 60 with another read (of unkonwn size) from the bank, then the percentage of shared kmer span is 75%. If a least a windows (of size \"windows_size\") contains at least kmer_threshold percent of position covered by shared kmers, the read couple is output."
    echo  "                 TRICK: with kmer_threshold<=0 a single kmer is sufficient in the linker mode to link two reads. "
    echo  "   -A index kmers present at least 'kmer_abundance_min' times in the bank AND in the queries."
    echo  "   -r (incompatible with SRC_counter), do not output precision about pair of similar reads. Only ids of reads from queries similar to at least one read from bank are output."
}


commet_like_option=""
bank_set=""
query_set=""
kmer_size=31
abundance_min=2
gamma=2
fingerprint_size=12
kmer_threshold=75
core_used=0
prefix="short_read_connector_res"
keep_low_complexity_option=""
remove=1
diskMode=0
abundanceMode=0 # 0: only bank, 1: bank and query"
countMode=0
windows_size=0

#######################################################################
#################### GET OPTIONS                #######################
#######################################################################
while getopts "hgb:q:p:k:a:s:t:f:G:w:dcrlA" opt; do
    case $opt in

    h)
        help
        exit
        ;;

    G)
        echo "use gamma: $OPTARG" >&2
        gamma=$OPTARG
        ;;

    d)
        echo "use disk mode">&2
        diskMode=1
        ;;

    r)
        echo "Output only ids of read shared (no complete links)">&2
        commet_like_option="-no_sharing_detail"
        ;;

    w)
        echo "use windows size: $OPTARG" >&2
        windows_size=$OPTARG
        ;;


    c)

        echo "use SRC_counter">&2
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
        echo "reuse precomputed solid kmers if exists"
        remove=0
        ;;
    A)
        echo "kmers must be present at least \'a\' times in the bank AND in the queries to be indexed"
        abundanceMode=1
        ;;

    l)
        echo "keep low complexity sequences"
        keep_low_complexity_option="-keep_low_complexity"
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
    exit 1
fi

# in case of linker (without the -r option that only checks for the presence of at least something in the bank), we cannot apply to multiple sequences
if [ $countMode -eq 0 ] && [ ! $commet_like_option ]; then
    # Check that the bank file is NOT a file of file. If its gzipped: consider its ok, else check that the first line contains something that is not a file
    gunzip -t ${bank_set} 2>/dev/null
    if [ $? -ne 0 ] # IN CASE THIS IS NOT A GZIPPED FILE
    then
        line=`head -n 1 ${bank_set} | cut -f 1 -d ' '`
        if [ -f $line ]; then
            echo "With SRC linker (without the -r option), the bank must not be composed of a file of file. It must be a fasta or fastq file gzipped or not"
            exit 1
        fi
    fi
fi




if [ -z "${query_set}" ]; then
    echo "You must provide a query read set (-q)"
    exit 1
fi



out_dsk=${prefix}"_solid_kmers_k"${kmer_size}".h5"
result_file=${prefix}".txt"


if [ $remove -eq 1 ]; then
    rm -f ${out_dsk}
fi

if [ $countMode -eq 1 ]; then
    if [ $commet_like_option ]; then
        echo "        ERROR: options -c and -r incompatibles"
        exit 1       
    fi
fi
# Count kmers using dsk if file absent
if [ ! -e ${out_dsk} ]; then
    if [ $abundanceMode -eq 1 ]; then
        ls ${bank_set} ${query_set} > FOF_FOR_DSK_REMOVE_ME_PLEASE.txt 
        cmd="${dsk_bin} -file FOF_FOR_DSK_REMOVE_ME_PLEASE.txt -kmer-size ${kmer_size} -abundance-min ${abundance_min} -out ${out_dsk} -nb-cores ${core_used} -solidity-kind all"
    else
        cmd="${dsk_bin} -file ${bank_set} -kmer-size ${kmer_size} -abundance-min ${abundance_min} -out ${out_dsk} -nb-cores ${core_used} -solidity-kind one"
    fi 
    echo ${cmd}
    ${cmd}
    if [ $? -ne 0 ]
    then
        echo "there was a problem with the kmer counting."
        exit 1
    fi
fi



# Compare read sets


# SRC_LINKER_RAM
if [ $diskMode -eq 0 ]; then
    if [ $countMode -eq 0 ]; then
        cmd="${BIN_DIR}/SRC_linker_ram"
        threshold_option="-kmer_threshold ${kmer_threshold}"
    else
        # SRC_COUNTER
        cmd="${BIN_DIR}/SRC_counter"
        threshold_option="-coverage_threshold ${kmer_threshold}"
    fi
else
    echo "Disk version not maintained anymore - sorry"
    exit 1;

    # SRC_LINKER_DISK
    #cmd="${BIN_DIR}/SRC_linker_disk"
fi


# Using a non-zero high density window: uncomment this line:
# zerod_option="-zero_density_windows_size 50 -zero_density_threshold 20"


# adding options
cmd="${cmd} -graph ${out_dsk}  -bank ${bank_set} -query ${query_set} -out ${result_file} ${threshold_option} -fingerprint_size ${fingerprint_size} -core ${core_used} -gamma ${gamma} ${commet_like_option} ${zerod_option} ${keep_low_complexity_option}"




# adding windows size option in the linker case
# if [ $countMode -eq 0 ]; then
    cmd="${cmd} -windows_size ${windows_size}"
# fi
echo ${cmd}
${cmd}
if [ $? -ne 0 ]
then
    echo "there was a problem with short read connector."
    exit 1
fi



# sort results
#sort -n ${unsorted_result_file} > ${result_file}


#rm -f ${unsorted_result_file}


echo "***********************************"
echo "Short read connector finished"
echo "results in:"
echo "   "${result_file}
echo "Contact: pierre.peterlongo@inria.fr"
echo "***********************************"
