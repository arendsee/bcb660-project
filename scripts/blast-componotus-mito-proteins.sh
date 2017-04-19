module load python/3.6.0
module load ncbi-blast

[[ -d smof ]] || git clone https://github.com/arendsee/smof

name='Camponotus atrox'
all_pro='data_mitochondria/mitochondrion*protein.faa'

my_pro=data_mitochondria/camponotus-mito.faa

cat $all_pro | smof/smof.py grep "$name" > $my_pro

tblastn                  \
    -task tblastn-fast   \
    -query $my_pro       \
    -db blastdb/bioblitz \
    -evalue 1e-2         \
    -num_threads 8       \
    -outfmt '7 qseqid sseqid qstart qstop sstart sstop bitscore evalue' > camponotus-blast.tab
