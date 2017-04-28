#!/usr/bin/env bash

usage (){
cat << EOF
USAGE
  scripts/make-histogram.sh <blast result> 
ARGUMENTS
  -n index of column with the species name (default=2)
  -b index of column with bitscore (default=8)
  -h show this help message
OUTPUT
  1. *-max-scores.tab - table mapping species to max bitscores
  2. *-plot.pdf - histogram of scores
DEPENDS
  Rscript
EXAMPLE
  scripts/make-histogram.sh nd6-blast-with-names.tab
EOF
    exit 0
}

# print help with no arguments
[[ $# -eq 0 ]] && usage

name_column=2
bitscore_column=8
while getopts "hn:b:" opt; do
    case $opt in
        n)
            name_column=$OPTARG ;;
        b)
            bitscore_column=$OPTARG ;;
        h)
            usage ;;
    esac 
done

input=$1
shift

out_table=${input/.tab/-max-scores.tab}
out_plot=${input/.tab/-plot.pdf}

if [[ ! -f $out_table ]]
then
    awk -v SCORE=$bitscore_column     \
        -v NAME=$name_column          \
    '
        BEGIN{ FS="\t" ; OFS="\t" }

        $SCORE > a[$NAME] { a[$NAME] = $SCORE }

        END { for (k in a){ print k, a[k] } }
    ' $input > $out_table
fi

Rscript - $out_table $out_plot << EOF
args <- commandArgs(TRUE)
d <- read.table(file=args[1])
pdf(args[2])
hist(d\$V2)
dev.off()
EOF
