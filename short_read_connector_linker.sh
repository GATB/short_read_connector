#!/bin/bash

version="1.2.0"



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
    echo "short_read_connector_linker.sh - Compare reads from two read sets (distinct or not)"
    echo "Version "$version

    echo "Usage: sh short_read_connector_linker.sh [index/query]"
}




function help_index {
    echo "short_read_connector_linker.sh index - Index a read file"
    echo "Usage: sh short_read_connector_linker.sh index -b read_file  -i dumped_index_name [OPTIONS]"
    echo  "   -b read_files for bank"
    echo  "     Example: -b data/c1.fasta.gz"

    echo  "   -i <string>. File of the index file to be created. Example \"my_index.dumped\""
    echo  "   -k <int> value (<32). Set the length of used kmers. Must fit the compiled value. Default=31"
    echo  "   -f <int> value. Fingerprint size. Size of the key associated to each indexed value, limiting false positives. Default=12"
    echo  "   -G <int> value. gamma value. MPHF expert users parameter - Default=2"
    echo  "   -a <int> kmer_abundance_min (kmer from bank seen less than this value both in the bank are not indexed). Default=2"
    echo  "   -l Keep low complexity regions (default false)"
    echo  "   -t <int> number of thread used. Default=0 (all)"
}


function help_query {
    
    echo "Usage: sh short_read_connector_linker.sh query -i indexed_file -q read_file_of_files [OPTIONS]"
    echo  "MANDATORY:"
    echo  " -i index created by short_read_connector_linker.sh index (read only)"
    echo  "   Example: -i my_index.dumped"
    echo  " -q read_file_of_files for query"
    echo  "   Example: -q data/fof.txt (with fof being a file of file descriptor)"

    echo  "OPTIONS:"
    echo  "   -w <int> window_size. See option -s. If the windows size is zero (default value), then the full read is considered. Default=0"
    echo  "   -s <int> kmer_threshold: Minimal percentage of shared kmer span for considering 2 reads as similar.  "
    echo  "            The kmer span is the number of bases from the read query covered by a kmer shared with the target read."
    echo  "            If a read of length 80 has a kmer-span of 60 with another read (of unkonwn size) from the bank, then the percentage of shared kmer span is 75%. If a least a windows (of size \"windows_size\") contains at least kmer_threshold percent of position covered by shared kmers, the read couple is output."
    echo  "                 TRICK: with kmer_threshold<=0 a single kmer is sufficient in the linker mode to link two reads. "
    echo  "   -l Keep low complexity regions (default false)"
    echo  "   -p <string> prefix. All out files will start with this prefix. Default=\"short_read_connector_res\""
    echo  "   -t <int> number of thread used. Default=0 (all)"



    # echo  "   -A index kmers present at least 'kmer_abundance_min' times in the bank AND in the queries."
    echo  "   -r do not output precision about pair of similar reads. Only ids of reads from queries similar to at least one read from bank are output."
}


if [ "$1" != "index" ] && [ "$1" != "query" ]; then 
    help
    exit 1
fi

feature=$1
shift

# #####
# INDEX
# #####
if [ "$feature" == "index" ]; then
    
    bank_set=""
    index_name=""
    kmer_size=31
    abundance_min=2
    gamma=2
    fingerprint_size=12
    core_used=0
    keep_low_complexity_option=""
    remove=1
    #######################################################################
    #################### GET INDEX OPTIONS          #######################
    #######################################################################
    while getopts ":hG:i:f:b:lk:a:t:" opt; do
        case $opt in

        h)
            help_index
            exit
            ;;

        G)
            echo "use gamma: $OPTARG" >&2
            gamma=$OPTARG
            ;;
        i)
            echo "create index $OPTARG" >&2
            index_name=$OPTARG
            ;;

        f)
            echo "use fingerprint size: $OPTARG" >&2
            fingerprint_size=$OPTARG
            ;;

        b)
            echo "use bank read set: $OPTARG" >&2
            bank_set=$OPTARG
            ;;
            

        l)
            echo "keep low complexity sequences"
            keep_low_complexity_option="-keep_low_complexity"
            ;;

        k)
            kmer_size=$OPTARG
            if [ $kmer_size -gt 31 ]; then
                echo "ERROR, kmer size must be <32. Please choose a lower k value"
                exit 1
            fi
            echo "use k=$OPTARG" >&2
            ;;


        a)
            echo "use abundance_min=$OPTARG" >&2
            abundance_min=$OPTARG
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
    #################### END GET INDEX OPTIONS      #######################
    #######################################################################

    if [ -z "${bank_set}" ]; then
        echo "You must provide a bank read set (-b)"
        help_index
        exit 1
    fi

    if [ -z "${index_name}" ]; then
        echo "You must provide an index name (-i)"
        help_index
        exit 1
    fi


    # in case of linker, we cannot apply to multiple sequences
    if [ ! $commet_like_option ]; then
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


    out_dsk="solid_kmers_k"${kmer_size}".h5"
# echo $abundanceMode
# exit
    # if [ "$abundanceMode" -eq 1 ]; then
    #     ls ${bank_set} ${query_set} > FOF_FOR_DSK_REMOVE_ME_PLEASE.txt 
    #     cmd="${dsk_bin} -file FOF_FOR_DSK_REMOVE_ME_PLEASE.txt -kmer-size ${kmer_size} -abundance-min ${abundance_min} -out ${out_dsk} -nb-cores ${core_used} -solidity-kind all"
    # else
        cmd="${dsk_bin} -file ${bank_set} -kmer-size ${kmer_size} -abundance-min ${abundance_min} -out ${out_dsk} -nb-cores ${core_used} -solidity-kind one"
    # fi 
    echo ${cmd}
    ${cmd}
    if [ $? -ne 0 ]
    then
        echo "there was a problem with the kmer counting."
        exit 1
    fi

    # create the index
    cmd="${BIN_DIR}/SRC_linker  -make_index -dumped_quasi_dict ${index_name} -graph ${out_dsk} -bank ${bank_set} -fingerprint_size ${fingerprint_size} -core ${core_used} -gamma ${gamma}  ${keep_low_complexity_option}"

    echo ${cmd}
    ${cmd}
    if [ $? -ne 0 ]
    then
        echo "there was a problem with the index creation."
        exit 1
    fi

    # rm -f ${out_dsk}

    echo "***********************************"
    echo "Short read connector indexation finished"
    echo "index in:"
    echo "   "${index_name}
    echo "Contact: pierre.peterlongo@inria.fr"
    echo "***********************************"

fi # END INDEXATION

# #####
# QUERY 
# #####
if [ "$feature" == "query" ]; then
    commet_like_option=""
    query_set=""
    kmer_size=31
    kmer_threshold=75
    core_used=0
    prefix="short_read_connector_res"
    keep_low_complexity_option=""
    windows_size=0

    #######################################################################
    #################### GET QUERY OPTIONS          #######################
    #######################################################################
    while getopts "hi:rw:q:p:ls:t:" opt; do
        case $opt in
        h)
            help_query
            exit
            ;;


        i)
            echo "Use index $OPTARG" >&2
            index_name=$OPTARG
            ;;

        q)
            echo "use query read set: $OPTARG" >&2
            query_set=$OPTARG
            ;;


        w)
            echo "use windows size: $OPTARG" >&2
            windows_size=$OPTARG
            ;;

        p)
            echo "use prefix=$OPTARG" >&2
            prefix=$OPTARG
            ;;

        l)
            echo "keep low complexity sequences"
            keep_low_complexity_option="-keep_low_complexity"
            ;;

        t)
            echo "use $OPTARG threads">&2
            core_used=$OPTARG
            ;;

        s)
            echo "use kmer_threshold=$OPTARG" >&2
            kmer_threshold=$OPTARG
            ;;

        r)
            echo "Output only ids of read shared (no complete links)">&2
            commet_like_option="-no_sharing_detail"
            ;;

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
    #################### END GET QUERY OPTION       #######################
    #######################################################################


    if [ -z "${index_name}" ]; then
        echo "You must provide the index name (-i)"
        help_query
        exit 1
    fi


    if [ -z "${query_set}" ]; then
        echo "You must provide a query read set (-q)"
        help_query
        exit 1
    fi



    
    result_file=${prefix}".txt"

   



    # Using a non-zero high density window: uncomment this line:
    # zerod_option="-zero_density_windows_size 50 -zero_density_threshold 20"


    # Creating the command
    cmd="${BIN_DIR}/SRC_linker -dumped_quasi_dict ${index_name} -query ${query_set} -out ${result_file} -kmer_threshold ${kmer_threshold}  -core ${core_used} ${commet_like_option} ${zerod_option} ${keep_low_complexity_option} -windows_size ${windows_size}"



    echo ${cmd}
    ${cmd}
    if [ $? -ne 0 ]
    then
        echo "there was a problem with short read connector."
        exit 1
    fi


    echo "***********************************"
    echo "Short read connector finished"
    echo "results in:"
    echo "   "${result_file}
    echo "Contact: pierre.peterlongo@inria.fr"
    echo "***********************************"

fi

