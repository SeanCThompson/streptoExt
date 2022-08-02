#!/bin/sh

# Extract last row x col of data file

runName=$(basename $1 .txt | awk -F '[_]' '{print$2}')

tail -n 90000 ../src/data_$runName.txt > ${runName}_lastframe.txt

echo $runName

# Find most common lineage

lineage=$(awk -F '[ ]' '{print $4}' ${runName}'_lastframe.txt' | uniq -c | sort -nr | head  -1 | awk '{ gsub(/^[ \t]+|[ \t]+$/, ""); print }' | awk -F '[ ]' '{print $1}')

echo $lineage

genome=$(grep " $lineage " ${runName}'_lastframe.txt' | head -1 | awk '{$1=$2=$3=$4=""; print $0}')

#genome='AAAAAAAAAFFFFFFFFFFFFFFFFFFFFFAFFFFFFFFFFFFFFFFFFFFFFFFFA 38082,38018,50212,1154,34086,38082,38018,34214,50212,50214,50342, 0.001 0.8, 0.9 0.9, 1 1'
#genome='AFAFFFFAFAFFFFFFFFFFFAFFFAFFFFFAFFFFFFFFFFFAFFAFAFFFFFFFFFAFF 38082,38018,50212,1154,34086,38082,38018,34214,50214,50214,50342, 0.046222 0.811903, 0.883708 0.618599, 1 1'
#genome='AAAAAAAAAAAFFFFFFFFFFFFFFFFFFFF 38082,38018,50212,1154,34086,38082,38018,34214,50212,50214,50342, 0.001 0.8, 0.9 0.9, 1 1'
genome='FFFFFFFFFFFAFFAFAFFFFFFFFAAFFFFFAFFFFFAFAFFFFAAAFA 47426,47939,47427,47939,47427,49006,15593,15465,48963,47939,48361,12365, 0.0104183 0.832070, 0.984682 0.737399, 1 1'


echo $genome



echo ${runName}_logging

../src/./strepto -name ${runName}_logging -v 1 -breakpoint_mut_type P -breakpoint_init 0 -breakprob 0 -breakpoint_inflow 0 -stressed_break 0.05 -ddrate 0 -maxtime 2500 -season 50001 -trans_prob_mut 0.0 -c 300 -r 300 -prob_mut_antibtype_perbit 0.0 -contact_break 1 -init_genome $genome 

#python3 correlate_ABorR_Stress.py ab_mut_log_${runName}_logging.txt

python3 correlate_ABorR_intime.py ab_mut_log_${runName}_logging.txt


