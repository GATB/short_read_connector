#!/bin/bash

sh ../short_read_connector.sh -b ../data/c1.fasta.gz -q ../data/c2.fasta.gz -p linker
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on linker"
  exit 1
fi

sort linker.txt > linker_sorted.txt
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on sorting linker"
  exit 1
fi


diff linker_sorted.txt ref_linker_sorted.txt
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on diff linker"
  exit 1
fi



sh ../short_read_connector.sh -b ../data/c1.fasta.gz -q ../data/c2.fasta.gz -d -p linker_disk
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on linker_disk"
  exit 1
fi

sort linker_disk.txt > linker_disk_sorted.txt
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on sorting linker_disk"
  exit 1
fi


diff linker_disk_sorted.txt ref_linker_sorted.txt
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on diff linker_disk"
  exit 1
fi



sh ../short_read_connector.sh -b ../data/c1.fasta.gz -q ../data/c2.fasta.gz -c -p counter
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on counter"
  exit 1
fi

sort counter.txt > counter_sorted.txt
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on sorting counter"
  exit 1
fi

diff counter_sorted.txt ref_counter_sorted.txt
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on diff counter"
  exit 1
fi

rm -f Erase_Me *.h5 counter.txt counter_sorted.txt linker.txt linker_disk.txt linker_disk_sorted.txt linker_sorted.txt

echo "*** Test: OK"

exit 0