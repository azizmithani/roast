# ROAST: a tool for reference-free optimization of supertranscriptome assemblies

Program: ROAST (Reference-free Optimization of Assembled SuperTranscriptomes)

Version: 1.0.0 (using bamtools 2.4.0 and boost c++ libraries)

## Dependencies:

ROAST is written in C++ and needs the following programs (provided as external tools with ROAST):

* CD-HIT-EST

* MINIMAP2

* HISAT2

* Cufflinks

* BLAT

* CAP3

* Picard

* supertranscript script 

Following programming languages and tools must be installed on the system and added in PATH (required)

* Java

* Python

* C++

* BLAST

* Samtools (version >= 1.9)

* BOOST API library

* Bamtools API library (github.com/pezmaster31/bamtools)

* Trinity (required if assembly is not provided)


Alternatively, create Conda environment and install all pre-requiste libraries and tools from the roast-env.sh file using following commands

	conda create --name roast-env python=3.8 --yes
	conda activate roast-env

Now run roast-env.sh from the terminal to install all the required tools and libraries.

**Bamtools and BOOST library Installation:**

Set library path

	export LD_LIBRARY_PATH="/path/to/bamtools/installation/dir/lib:/path/to/boost/installation/dir/lib"
	
Or if you have used conda environment for installation, set the library path as:

	export LD_LIBRARY_PATH="/path/to/conda/envs/roast-env/lib"


## Commands to Compile:
To compile ROAST in the working directory, copy and paste the following bash commands in a text file and save as a bash script. Alternatively, run the roast-compile.sh script from the working directory of ROAST. Before running this bash script, please change the paths for boost and bamtools installation directories with yours.

	====
	#!/bin/bash
	bamtools=/path/to/my/bamtools
	boost=/path/to/my/boost
	
	 ====
#!/bin/bash
bamtools=/path/to/my/bamtools
boost=/path/to/my/boost

 unzip external_tools/cufflinks.zip -d external_tools/

 g++ -I $bamtools/include/bamtools -I $boost/include/ -L $boost/lib/ -L $bamtools/lib/ ROAST_extendContigs.cpp -o ROAST_extendContigs -lbamtools  -lboost_filesystem -lboost_regex -lz
 g++ -I $bamtools/include/bamtools -I $boost/include/ -L $boost/lib/ -L $bamtools/lib/ ROAST_extendContigs_SCs.cpp -o ROAST_extendContigs_SCs  -lbamtools -lboost_filesystem -lboost_regex -lz
 g++ -I $bamtools/include/bamtools -I $boost/include/ -L $boost/lib/ -L $bamtools/lib/ ROAST_mergeContigs_SCs.cpp -o ROAST_mergeContigs_SCs -lbamtools -lboost_filesystem -lboost_regex -lz
 g++ -o roast -I $bamtools/include/bamtools -I $boost/include/ -L $boost/lib/ -L $bamtools/lib/ main.cpp mis_assembly_chimera.cpp global.cpp filterSamToFastq.cpp bySoftclip.cpp byRI.cpp alignment.cpp -lbamtools -lboost_filesystem -lboost_regex -lz

====
	
	 g++ -I $bamtools/include/bamtools -I $boost/include/ -L $boost/lib/ -L $bamtools/lib/ ROAST_extendContigs.cpp -o ROAST_extendContigs -lbamtools  -lboost_filesystem -lboost_regex -lz
	 g++ -I $bamtools/include/bamtools -I $boost/include/ -L $boost/lib/ -L $bamtools/lib/ ROAST_extendContigs_SCs.cpp -o ROAST_extendContigs_SCs  -lbamtools -lboost_filesystem -lboost_regex -lz
	 g++ -I $bamtools/include/bamtools -I $boost/include/ -L $boost/lib/ -L $bamtools/lib/ ROAST_mergeContigs_SCs.cpp -o ROAST_mergeContigs_SCs -lbamtools -lboost_filesystem -lboost_regex -lz
	 g++ -o roast -I $bamtools/include/bamtools -I $boost/include/ -L $boost/lib/ -L $bamtools/lib/ main.cpp mis_assembly_chimera.cpp global.cpp filterSamToFastq.cpp bySoftclip.cpp byRI.cpp alignment.cpp -lbamtools -lboost_filesystem -lboost_regex -lz
	
	====


ROAST can be run from installation/directory/ROAST or by setting the environmental variable ROAST to point to this, which will make it easy to access both ROAST as well as external tools that come bundled with ROAST.

 	export ROAST=/path/to/ROAST/installation/dir

 	$ROAST/roast --help (for help)


## Usage:

	roast --fastq_1 fastq1_filename --fastq_2 fastq2_filename <Reference type> --threads INT [Parameters]

**Reference type:**

	--supertranscript_assembly <STR> 	Provide SuperTranscript fasta for supertranscriptome assembly improvement

	--trinity_assembly <STR>		Provide transcript fasta for raw Trinity assembly improvement

	--generate_assembly 			Generate de novo supertranscriptome assembly for improvement. Set and export the environmental variable TRINITY_HOME to point to the Trinity installation folder 

**Parameters:**
   
	--output_dir <STR>		 	Path for output directories, default folder of input reference sequence
   
	--inner_itr <INT> 			Number of inner iterations threshold, default 30
   
	--outer_itr <INT> 			Number of outer Iterations threshold, default 100

	--improvment_TH <INT> 		 	Keep improving until number of improved contigs meet threshold, default 1

	--min_extended_contigs <INT> 	 	Keep extending until number of extended contigs meet threshold, default 1

	--max_memory_TRINITY <INT> 		Maximum memory for TRINITY, default 20
	
	--threads <INT> 			Number of threads, default 8
   
	--threadsForSamSort <INT> 		Number of threads for samtools sort, default 2
	
	--memForSamSort <INT>			Memory for samtools sort, default 768M
   
	--complete_cleanup 			Complete cleanup of intermediate files, default partial cleanup
	
	--nochange_header 			Don't change header of final improved assembly file, default change header
	
	--cdHitEST <INT> 			0 for no CD-HIT-EST, 1 for CD-HIT-EST only in the start, 2 for CD-HIT-EST after every iteration, default 1
	
	--mapping_quality <INT>			Minimum mapping quality to filter Bam file, default 20
	
	--sc_support_TH <FLOAT>			Softclips support for a position for Incomplete and fragmented contigs and false chimera process, default 0.75
	
	--min_SCs_reads <INT> 		  	Minimum number of softclips reads for extension of Incomplete contigs, default 3
	
	--sc_cons_len_TH <INT> 		 	Minimum length of consensus sequence generated by softclips to BLAST for fragmented contigs Identification, default 10
	
	--discard_contig_corner_len_TH <INT>  	Minimum number of bases allowed to discard from the corner of contigs while finding overlapped edges, default 10
	
	--max_allowed_gaps <INT> 		Maximum gaps allowed to find overlaps between the corners of the two contigs, default 0
	
	--edge_boundary <INT>			Terminal region for mapped reads to be considered for un-mapped and distantly mapped reads, default 2x of read length
	
	--min_unmapped_reads <INT> 		Minimum number of unmapped reads to generate CAP3 assembly for incomplete contigs extension, default 5
	
	--min_distMapped_reads <INT> 	 	Minimum number of distantly mapped reads to consider as read island for merging of fragmented contigs, default 3
	
	--contig_boundary	 		Check overlap between two fragmented contigs within 5% read length of the contig boundary, default read island boundary
	
	--min_allowed_unmapped_ext <INT> 	Minimum number of bases to be considered for valid extension using unmapped reads, default 50% of read length
	
	--one_side_allowed_SCs_RI <INT> 	Maximum % of softclips of the total read length allowed at one side of read to consider it for Read Island, default 25
	
	--each_side_allowed_SCs_RI <INT> 	Maximum % of softclips of the total read length allowed at both sides of read to consider it for Read Island, default 12
	
	--sc_start_pos_from_terminus <INT> 	For terminal softclips extraction, define starting position to consider, default 25
	
	--min_overlap_TH <INT> 		 	Minimum length of overlapped sequence to consider for merging fragmented contigs and CAP3 assemblies, default 20 bases
	
	--win_size <INT> 			Window size to detect gradual coverage change, default 2X read_length
	
	--win_diff_TH <INT> 			Maximum average coverage change threshold between two consecutive windows, default off
	
	--coverage_drop_TH <FLOAT> 		Maximum coverage ratio between two consecutive positions to process for false chimera identification process, default 0.2
	
	--st_end_boundary <INT> 		Consider coverage change within start and end boundary of the contig, default one and half of read length 
	
	--blast_score_TH <INT> 		 	BLAST hit identity and coverage score to find overlap between fragmented contigs and identify mis-assembly/false chimera, default 90" << endl;

	
	--ignore_short_seq <INT> 		Minimum length of the contig to remove from the assembly, default 200
	
	--insert_Ns <INT> 			Number of Ns to insert between contig and its CAP3 assembly in the absence of overlap, default 5
	
**Terminate ROAST process:**

To stop ROAST properly at the end of current iteration, place an empty file named 'stop.txt' in the folder 'intermediate_Improved_assemblies'.

   
## Output:

Improved fasta files from all iterations named as Name_SuperIterationNumber-MiniIterationNumber.fasta e.g. Name_1-1.fasta for first improvement in intermediate_improved_assemblies folder.

Final improved file named final_improved.fasta, summary file and removed sequences can be found in final_improved assembly folder.

The log data for each improvement at each step is present in log folder with file name after each super and mini iterations. Time log shows time taken by each step.

Run statistics are displayed on screen.

## Test data:

Test data is available here:
https://zenodo.org/record/8192067

Arabidopsis test dataset consists of paired-end fast files, initial assembly in fasta format, improved assembly by ROAST along with log files and summary output, and TransRate scoring matrices for both initial and improved assembly.

## Citing ROAST:

Please cite as

Shabbir, M., Mithani, A. Roast: a tool for reference-free optimization of supertranscriptome assemblies. BMC Bioinformatics 25, 2 (2024). https://doi.org/10.1186/s12859-023-05614-4
