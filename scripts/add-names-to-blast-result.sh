#!/usr/bin/bash

set -u

# USAGE:
#   # assuming you are running this from the bcb660-project folder:
#   scripts/autoslurm.sh scripts/add-names-to-blast-result.sh nd6-blast.tab

# This should be the name of the output file from BLAST, e.g. nd6-blast.tab
# 1. The file may have comments inside
# 2. The file is expected to have the extension `.tab`
input_blast_result=$1

# --- You shouldn't need to change this
# The protein fasta files are used just to extract the map between reference
# and species name. The fasta headers are assumed to have the format:
# >REFERENCE ... [ SPECIES_NAME ]
protein_reference=data_mitochondria/mitochond*protein.faa

# This is the name of the output file, the parameter expansion below replaces
# the extension on the input.
output_blast_result=${input_blast_result/.tab/-with-names.tab}

# Extract a reference to name map from the protein sequence file
# Join the map table with the blast result table
#  1. This adds the scientific species names as the second column in the output file
#  2. Replaces space in the species name with underscores
#  3. This will remove all comments from the blast result file 
join -t $'\t' \
  <( cat $protein_reference                   |
      sed -nr 's/>([^ ]+).*\[(.*)\]/\1\t\2/p' |
      tr ' ' '_' | sort ) \
  <( grep -Pv '^#' $input_blast_result | sort ) > $output_blast_result
