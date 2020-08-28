#!/bin/bash

##########################################################
## TEST RAM LINKER
##########################################################

# RUN SRC
# Create the index
(bash ../short_read_connector_linker.sh index -b ../data/c1.fasta.gz -i index_linker.dumped)> log_linker 2> log_linker_err
# Permorf queries
(bash  ../short_read_connector_linker.sh query -i index_linker.dumped -q fof.txt -p linker)>> log_linker 2>> log_linker_err
# (bash ../short_read_connector_linker.sh -b ../data/c1.fasta.gz -q fof.txt -p linker) > log_linker 2> log_linker_err
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
## TEST RAM LINKER WITHOUT LINKS
##########################################################

# RUN SRC

# Create the index
(bash ../short_read_connector_linker.sh index -b ../data/c1.fasta.gz -i index_linker.dumped)> log_linker 2> log_linker_err
# Permorf queries
(bash  ../short_read_connector_linker.sh query -i index_linker.dumped -q fof.txt -p linker_no_link -r)>> log_linker 2>> log_linker_err


# (bash ../short_read_connector.sh -b ../data/c1.fasta.gz -q fof.txt -p linker_no_link -r) > log_linker 2> log_linker_err
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on linker without link output"
  exit 1
fi

cat linker_no_link.txt | wc |  sed 's/^[ \t]*//;s/[ \t]*$//' > wc_linker_no_link # sed for removing spaces at the beggin and at the end of the line

# CHECK EQUALITY
diff wc_linker_no_link ref_wc_linker_no_link
if [ $? -ne 0 ] ; then
  echo "*** Test: FAILURE on diff linker without link output"
  exit 1
fi

echo "*** DIFF LINKER WITHOUT LINKS OK ***"


##########################################################
## TEST COUNTER
##########################################################

# RUN SRC

# Create the index
(bash ../short_read_connector_counter.sh index -b ../data/c1.fasta.gz -i index_counter.dumped -t 1)> log_counter 2> log_counter_err
# Permorf queries
(bash  ../short_read_connector_counter.sh query -i index_counter.dumped -q fof.txt -p counter -t 1)>> log_counter 2>> log_counter_err

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

rm -f Erase_Me *.h5 counter.txt wc_* linker.txt linker_no_link.txt  index_counter.dumped index_linker.dumped

echo "*** Test: OK ***"

exit 0
