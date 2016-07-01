#!/bin/bash

##########################################################
## TEST RAM LINKER
##########################################################

# RUN SRC
(bash ../short_read_connector.sh -b ../data/c1.fasta.gz -q fof.txt -p linker) > log_linker 2> log_linker_err
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on linker"
  exit 1
fi

cat linker.txt | wc |  sed 's/^[ \t]*//;s/[ \t]*$//' > wc_linker # sed for removing spaces at the beggin and at the end of the line

# CHECK EQUALITY
diff wc_linker ref_wc_linker
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on diff linker"
  exit 1
fi

echo "*** DIFF LINKER OK ***"

##########################################################
## TEST DISK LINKER --> UNAVAILABLE FOR NOW
##########################################################

# # RUN SRC
# bash ../short_read_connector.sh -b ../data/c1.fasta.gz -q ../data/c2.fasta.gz -d -p linker_disk
# if [ $? -ne 0 ] ; then
#   echo "*** Test: FAILURE on linker_disk"
#   exit 1
# fi
#
# # UNIFORMIZE RESULTS
# python ../scripts/uniformize_SRC_linker_output.py linker_disk.txt > linker_disk_uniform.txt
# if [ $? -ne 0 ] ; then
#   echo "*** Test: FAILURE on uniformisator"
#   exit 1
# fi
#
# # SORT RESULTS
# sort -n linker_disk_uniform.txt > linker_disk_sorted.txt
# if [ $? -ne 0 ] ; then
#   echo "*** Test: FAILURE on sorting linker_disk"
#   exit 1
# fi
#
# # CHECK EQUALITY
# diff linker_disk_sorted.txt ref_linker_sorted.txt
# if [ $? -ne 0 ] ; then
#   echo "*** Test: FAILURE on diff linker_disk"
#   exit 1
# fi


##########################################################
## TEST COUNTER
##########################################################

# RUN SRC
bash ../short_read_connector.sh -b ../data/c1.fasta.gz -q fof.txt -c -p counter -t 1 > log_counter 2>log_counter_err
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on counter ***"
  exit 1
fi

cat counter.txt | wc |  sed 's/^[ \t]*//;s/[ \t]*$//' > wc_counter  # sed for removing spaces at the beggin and at the end of the line

# CHECK EQUALITY
diff wc_counter ref_wc_counter
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on diff counter"
  exit 1
fi

echo "*** DIFF COUNTER OK ***"

##########################################################
## CLEAN TEMP FILES
##########################################################

rm -f Erase_Me *.h5 counter.txt wc_* linker.txt log*

echo "*** Test: OK ***"

exit 0
