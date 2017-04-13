mitodir=data_mitochondria
mitorep=mitochondria-report.txt

species_list=species-list.txt

>$mitorep

count_them () {
    echo -n "Number of represented mitochondrial genomes: "
    zcat $mitdir/*protein*faa.gz |
        sed 's/.*\[\(.*\)\]$/\1/' |
        uniq | sort -u | tee $species_list | wc -l
}

count_genes () {
    echo "Number of genomes with n proteins:"
    echo -e "count\tn_proteins"
    zcat $mitdir/*protein*faa.gz |
        sed 's/.*\[\(.*\)\]$/\1/' |
        sort | uniq -c | tr ' ' "\t" | sort -rnk2
}

list_irregular_species_names () {
    echo "Species with irregular names (i.e. not '[A-Z][a-z]+ [a-z]+'):"
    grep -P '^[A-Z][a-z]+ [a-z]+$' $species_list | sed 's/^/  /'
}

count_them >> $mitorep
count_genes >> $mitorep
list_irregular_species_names >> $mitorep
