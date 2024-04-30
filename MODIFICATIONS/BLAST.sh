#!/bin/bash

# Load BLAST+ module, if necessary
source /apps/profiles/modules_asax.sh.dyn
module load blast+/2.15.0

# Define the database and query file paths
database="refseqgene_db"
query_file="/scratch/aubclsc0332/Lizard_alt_ref_Phrynosoma/LizardRefGenome/genes.fasta"

# Define the output file path
output_file="/scratch/aubclsc0332/Lizard_alt_ref_Phrynosoma/LizardRefGenome/blast_results.txt"

# Run BLASTn
blastn -query $query_file \
       -db $database \
       -out $output_file \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
       -max_target_seqs 1 \
       -evalue 0.001

# Output completion message
echo "BLAST search has completed. Results are in: $output_file"
