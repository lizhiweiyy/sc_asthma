#!/bin/bash
export PATH=~/software/cellranger-4.0.0:$PATH
for i in {"B6","E4","E5"}; do
cd ~/path/${i}/cleandata/
echo $i
pathfq=~/path/${i}/cleandata/
cellranger count --id=${i} \
	--sample=${i} \
	--transcriptome=Ref/refdata-gex-mm10-2020-A\
	--fastqs=$pathfq
done