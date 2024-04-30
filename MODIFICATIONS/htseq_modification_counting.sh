#!/bin/bash
# Load necessary modules
source /apps/profiles/modules_asax.sh.dyn
module load htseq
MyID=aubclsc0334
CountHT=/scratch/$MyID/PracticeRNAseq_Full_Script/Counts_HTSeq_modified
MAPD=/scratch/$MyID/PracticeRNAseq_Full_Script/Map_HiSat2_Old
REFD=/scratch/$MyID/PracticeRNAseq_Full_Script/LizardRefGenome
RESULTSD=/home/$MyID/PracticeRNAseq_Full_Script/Counts_H_S_2024

mkdir $CountHT
cd $CountHT
git clone https://github.com/htseq/htseq
cp $CountHT/htseq/HTSeq/scripts/count.py .

# List of variables
variables=("SRR6825925"  "SRR6825928"  "SRR6825932"  "SRR6825934"  "SRR6825938"  "SRR6825940" "SRR6825926"  "SRR6825930"  "SRR6825933"  "SRR6825935"  "SRR6825939" "SRR6825941")
# Iterate over each variable
for var in "${variables[@]}"
do
   # Create directory
   mkdir "$var"
   # Run HTSeq count
python -m HTSeq.scripts.count -f bam --max-reads-in-buffer=35000000 -a 20 -t exon -i gene_id --additional-attr=gene_name --additional-attr=exon_number -c $CountHT/"$var"/"$var".csv -o $CountHT/"$var"/"$var".gtf $MAPD/"$var"_sorted.bam $REFD/SceUnd1.0_top24.gtf > $CountHT/"$var"/"$var".tab
done
