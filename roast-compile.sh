#!/bin/bash
bamtools=~/miniconda2/envs/roast-env2
boost=~/miniconda2/envs/roast-env2

unzip external_tools/cufflinks

g++ -I $bamtools/include/bamtools -I $boost/include/ -L $boost/lib/ -L $bamtools/lib/ ROAST_extendContigs.cpp -o ROAST_extendContigs -lbamtools -lboost_filesystem -lboost_regex -lz

g++ -I $bamtools/include/bamtools -I $boost/include/ -L $boost/lib/ -L $bamtools/lib/ ROAST_extendContigs_SCs.cpp -o ROAST_extendContigs_SCs -lbamtools -lboost_filesystem -lboost_regex -lz

g++ -I $bamtools/include/bamtools -I $boost/include/ -L $boost/lib/ -L $bamtools/lib/ ROAST_mergeContigs_SCs.cpp -o ROAST_mergeContigs_SCs -	lbamtools -lboost_filesystem -lboost_regex -lz

g++ -o roast -I $bamtools/include/bamtools -I $boost/include/ -L $boost/lib/ -L $bamtools/lib/ main.cpp mis_assembly_chimera.cpp global.cpp filterSamToFastq.cpp bySoftclip.cpp byRI.cpp alignment.cpp -lbamtools -lboost_filesystem -lboost_regex -lz
