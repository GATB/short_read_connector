#!/bin/bash

##########################################################
## TEST RAM LINKER
##########################################################

# RUN SRC
bash ../short_read_connector.sh -b ../data/c1.fasta.gz -q ../data/c2.fasta.gz -p linker
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on linker"
  exit 1
fi

# UNIFORMIZE RESULTS
python ../scripts/uniformize_SRC_linker_output.py linker.txt > linker_uniform.txt
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on uniformisator"
  exit 1
fi

# SORT RESULTS
sort -n linker_uniform.txt > linker_sorted.txt
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on sorting linker"
  exit 1
fi

# CHECK EQUALITY
diff linker_sorted.txt ref_linker_sorted.txt
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on diff linker"
  exit 1
fi



##########################################################
## TEST DISK LINKER
##########################################################

# RUN SRC
bash ../short_read_connector.sh -b ../data/c1.fasta.gz -q ../data/c2.fasta.gz -d -p linker_disk
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on linker_disk"
  exit 1
fi

# UNIFORMIZE RESULTS
python ../scripts/uniformize_SRC_linker_output.py linker_disk.txt > linker_disk_uniform.txt
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on uniformisator"
  exit 1
fi

# SORT RESULTS
sort -n linker_disk_uniform.txt > linker_disk_sorted.txt
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on sorting linker_disk"
  exit 1
fi

# CHECK EQUALITY
diff linker_disk_sorted.txt ref_linker_sorted.txt
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on diff linker_disk"
  exit 1
fi


##########################################################
## TEST COUNTER
##########################################################

# RUN SRC
bash ../short_read_connector.sh -b ../data/c1.fasta.gz -q ../data/c2.fasta.gz -c -p counter
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on counter"
  exit 1
fi

# UNIFORMIZE RESULTS
grep -v "#" counter.txt > counter_uniform.txt

# SORT RESULTS
sort -n counter_uniform.txt > counter_sorted.txt
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on sorting counter"
  exit 1
fi

# CHECK EQUALITY
diff counter_sorted.txt ref_counter_sorted.txt
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on diff counter"
  exit 1
fi


##########################################################
## CLEAN TEMP FILES
##########################################################

rm -f Erase_Me *.h5 counter.txt counter_sorted.txt linker.txt linker_disk.txt linker_disk_sorted.txt linker_sorted.txt linker_uniform.txt linker_disk_uniform.txt counter_uniform.txt

echo "*** Test: OK"

exit 0
