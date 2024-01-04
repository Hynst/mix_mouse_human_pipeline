#!/bin/bash

RUN=$1
MYELO="/mnt/shared/MedGen/Myelo/runs"

mkdir -p ${MYELO}/$RUN/variants
mkdir -p ${MYELO}/$RUN/coverage

for csv in `ls /mnt/shared/S3bibs/data/mculen/$RUN/samples/*/variants/mutect2/*csv`
do

echo $csv
cp $csv ${MYELO}/$RUN/variants

done

for txt in `ls /mnt/shared/S3bibs/data/mculen/$RUN/samples/*/coverage/*perexon*`
do

echo $txt
cp $txt ${MYELO}/$RUN/coverage

done
