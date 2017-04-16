# Retrieve mitochondrial genomes
# If they are already downloaded, do nothing

outdir=data_mitochondria

if [[ ! -d $outdir ]]
then
    mkdir $outdir
    wget \
        -P $outdir -nd -r -A gz \
        ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion
    cd $outdir
    gunzip *gz
fi
