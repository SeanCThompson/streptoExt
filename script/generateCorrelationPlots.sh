#!/bin/sh

# Extract last row x col of data file

runName=$(basename $1 .txt | awk -F '[_]' '{print$2}')

tail -n 90000 data_$runName.txt > ${runName}_lastframe.txt

echo $runName

# Find most common lineage

lineage=$(awk -F '[ ]' '{print $4}' ${runName}'_lastframe.txt' | uniq -c | sort -nr | head  -1 | awk '{ gsub(/^[ \t]+|[ \t]+$/, ""); print }' | awk -F '[ ]' '{print $1}')

echo $lineage

genome=$(grep " $lineage " ${runName}'_lastframe.txt' | head -1 | awk '{$1=$2=$3=$4=""; print $0}')

echo $genome

echo ${runName}_logging

../src/./strepto -name ${runName}_logging -v 1 -init_genome $genome -ddrate 0 -maxtime 2000

python3 correlate_F_and_ABorR.py ab_mut_log_${runName}_logging.txt data_${runName}.txt 



