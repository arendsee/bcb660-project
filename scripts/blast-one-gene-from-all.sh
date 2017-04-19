module load python/3.6.0
module load ncbi-blast

[[ -d smof ]] || git clone https://github.com/arendsee/smof

protein_name="NADH dehydrogenase subunit 6"
all_pro='data_mitochondria/mitochondrion*protein.faa'

my_pro=data_mitochondria/nd6-mito.faa

cat $all_pro | smof/smof.py grep "$protein_name" > $my_pro

tblastn                  \
    -task tblastn-fast   \
    -query $my_pro       \
    -db blastdb/bioblitz \
    -evalue 1e-2         \
    -num_threads 8       \
    -max_target_seqs 1e8 \
    -outfmt '6 qseqid sseqid qlen slen qstart qstop sstart sstop bitscore evalue' > nd6-blast.tab

# ref|YP_009231648.1| NADH dehydrogenase subunit 2 (mitochondrion) [Camponotus atrox]
# ref|YP_009231649.1| cytochrome c oxidase subunit I (mitochondrion) [Camponotus atrox]
# ref|YP_009231650.1| cytochrome c oxidase subunit II (mitochondrion) [Camponotus atrox]
# ref|YP_009231651.1| ATP synthase F0 subunit 8 (mitochondrion) [Camponotus atrox] 0
# ref|YP_009231652.1| ATP synthase F0 subunit 6 (mitochondrion) [Camponotus atrox]
# ref|YP_009231653.1| cytochrome c oxidase subunit III (mitochondrion) [Camponotus atrox]
# ref|YP_009231654.1| NADH dehydrogenase subunit 3 (mitochondrion) [Camponotus atrox]
# ref|YP_009231655.1| NADH dehydrogenase subunit 5 (mitochondrion) [Camponotus atrox]
# ref|YP_009231656.1| NADH dehydrogenase subunit 4 (mitochondrion) [Camponotus atrox]
# ref|YP_009231657.1| NADH dehydrogenase subunit 4L (mitochondrion) [Camponotus atrox]
# ref|YP_009231658.1| NADH dehydrogenase subunit 6 (mitochondrion) [Camponotus atrox] 33
# ref|YP_009231659.1| cytochrome b (mitochondrion) [Camponotus atrox]
# ref|YP_009231660.1| NADH dehydrogenase subunit 1 (mitochondrion) [Camponotus atrox]
