#!/usr/bin/env bash

git clone https://github.com/arendsee/litade
git clone https://github.com/arendsee/smof

my_clade=6656 # taxon id for phylum Arthropoda
my_taxids=arthropoda-taxids.tab
my_scinames=arthropoda-scinames.tab
mito_faa=data_mitochondria/mitochondrion.*.protein.faa
my_faa=arthropoda-mito.faa

# Find taxon ids of all species in Arthropoda
perl litade/litade.pl $my_clade | cut -f2 > $my_taxids

# Get the taxonomy dump
wget -O taxdmp.zip ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
unzip -p taxdmp.zip names.dmp > names.dmp
rm taxdmp.zip

# Build a map from Taxon Id to Scientific name
sed 's/[\t ]*|[\t ]*/\t/g' names.dmp |
    awk '
        BEGIN{FS="\t"; OFS="\t"}
        /scientific name/ {print $1,$2}
    ' > taxid2sciname.tab

# Get the scientific name for every species in Arthropoda
join <(sort $my_taxids) <(sort taxid2sciname.tab) -t $'\t' | cut -f2 > $my_scinames

# Extract the protein sequence of all arthropods
smof/smof.py grep -w '\[([^\]]+)\]' -f $my_scinames $mito_faa > arthropod-mito.faa
