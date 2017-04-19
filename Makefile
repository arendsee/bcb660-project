.PHONY: run
run:
	[[ -d blastdb ]] || scripts/autoslurm.sh "--time=24:00:00" scripts/build_blast_db.pbs
	[[ -d data_mitochondria ]] || scripts/autoslurm.sh scripts/retrieve_mitochondria.sh
	[[ -f data_mitochondria/6656-mito.faa ]] || scripts/autoslurm.sh scripts/extract-arthropoda-mitochondria.sh

.PHONY: clean
clean:
	rm -rf blastdb

.PHONY: archive
archive:
	[[ -d archive ]] && mv archive archive-${RANDOM}
	mkdir archive
	[[ -d blastdb ]] && mv blastdb archive
