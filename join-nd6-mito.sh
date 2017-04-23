#!/usr/bin/bash

join <( cat data_mitochondria/mitochond*protein.faa           |
        sed -nr 's/>(ref.YP_003058230.1.).*\[(.*)\]/\1\t\2/p' |
        tr ' ' '\t' | sort )
     <( sort nd6-blast.tab ) > nd6-blast-with-names.tab
