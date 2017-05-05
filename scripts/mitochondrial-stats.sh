set -u

mitodir=data_mitochondria

statdir=mito-stats

[ -d $statdir ] || mkdir $statdir

species_list=$statdir/species-list.txt
gene_counts=$statdir/gene-counts.tab
irregular_names=$statdir/irregular-names.tab
mitorep=$statdir/mitochondria-report.txt

mitoprot=$mitodir/*protein*faa

list-species () {
    sed -nr 's/^>.*\[(.*)\]$/\1/p' $mitoprot |
        uniq | sort -u
}

count-genes () {
    echo -e "count\tn_proteins"
    sed -nr 's/^>.*\[(.*)\]$/\1/p' $mitoprot |
        sort | uniq -c |
        awk '{print $1}' |
        sort -n | uniq -c | sort -rn |
        sed 's/^ *//; s/  */\t/'
}

list-irregular-species-names () {
    grep -vP '^[A-Z][a-z]+ [a-z]+$' $species_list | sed 's/^/  /'
}

make-report () {
    echo -n "Number of represented mitochondrial genomes: "
    wc -l $species_list | sed 's/ .*//'

    echo
    echo "First 10 most common mitochondrial gene counts"
    head -11 $gene_counts

    echo
    echo "Species with irregular names (i.e. not '[A-Z][a-z]+ [a-z]+'):"
    wc -l $irregular_names | sed 's/ .*//'
}

[ -f $species_list    ] || list-species                 > $species_list
[ -f $irregular_names ] || list-irregular-species-names > $irregular_names
[ -f $gene_counts     ] || count-genes                  > $gene_counts
[ -f $mitorep         ] || make-report                  > $mitorep
