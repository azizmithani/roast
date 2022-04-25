# roast
ROAST: A tool for Reference-free Optimisation of Assembled SuperTranscriptomes
DEPENDENCIES:

ROAST is written in C++ and needs the following programs (provided as external tools with ROAST):

CD-HIT-EST
MINIMAP2
HISAT2
Cufflinks
BLAT
CAP3
Picard

Following programming languages and tools must have installed on the user system and add in PATH (required)

JAVA
python
c++
BLAST
Samtools (version >= 1.9)
BOOST API library
Bamtools API library (github.com/pezmaster31/bamtools)


Follwoing tools Install and add in PATH (if required) 
Trinity


INSTALLATION:

Set library path

	% export LD_LIBRARY_PATH="pathTo/installed/bamtools-2.4.0/lib:pathTo/installed/boost/stage/lib"

Command To Compile:

1.	g++ -I pathTo/installed/bamtools-2.4.0/src/ -I pathTo/installed/boost/INCLUDE/dir/ -L pathTo/boost_1_70_0/lib/ -L pathTo/bamtools-2.4.0/lib/ ROAST_extendContigs.cpp -o ROAST_extendContigs  -lbamtools -lboost_filesystem  -lboost_regex -lz

2.	g++ -I pathTo/installed/bamtools-2.4.0/src/ -I pathTo/installed/boost/INCLUDE/dir/ -L pathTo/boost_1_70_0/lib/ -L pathTo/bamtools-2.4.0/lib/  ROAST_extendContigs_SCs.cpp -o ROAST_extendContigs_SCs.cpp -o ROAST_extendContigs_SCs  -lbamtools -lboost_filesystem  -lboost_regex -lz

3.	g++ -I pathTo/installed/bamtools-2.4.0/src/ -I pathTo/installed/boost/INCLUDE/dir/ -L pathTo/boost_1_70_0/lib/ -L pathTo/bamtools-2.4.0/lib/  ROAST_extendContigs_SCs.cpp -o ROAST_mergeContigs_SCs  -lbamtools -lboost_filesystem  -lboost_regex -lz

4.	g++ -I pathTo/installed/bamtools-2.4.0/src/ -I pathTo/installed/boost/INCLUDE/dir/ -L pathTo/boost_1_70_0/lib/ -L pathTo/bamtools-2.4.0/lib/ main.cpp mis_assembly_chimera.cpp global.cpp filterSamToFastq.cpp bySoftclip.cpp byRI.cpp alignment.cpp utils.h -o roast -lbamtools -lboost_filesystem  -lboost_regex -lz


You can Run ROAST from installation/directory/ROAST or set the environmental variable ROAST to point to this, which will make it easy to access both ROAST as well as external tools that come bundled with ROAST.

 	% export ROAST=/path/to/ROAST/installation/directory

 	% $ROAST/roast --help (for help)


USAGE:

Program: ROAST (Reference free Optimization of Assembled Supertranscriptomes)
Version: 1.0.0 (using bamtools 2.4.0 and boost c++ libraries)


	roast --fastq_1 fastq1_filename --fastq_2 fastq2_filename <Reference type> --threads INT [Parameters]

Reference type:

	--supertranscript_assembly 		Provide SuperTranscript fasta for SuperTranscript assembly improvement

	--trinity_assembly 			Provide transcript fasta for Raw Trinity assembly improvement

	--generate_assembly 			Set and export the environmental variable TRINITY_HOME to point Trinity installation folder to generate De novo Transcriptome assembly and improve output

Parameters:
   
	--output_dir 			 	Path for output directories, default folder of input reference sequence
   
	--inner_itr <INT> 			Number of Inner iterations threshold, default 30
   
	--outer_itr <INT> 			Number of Outer Iterations threshold, default 100

	--improvment_TH <INT> 		 	Keep improving until number of improved contigs meet threshold, default 1

	--min_extended_contigs <INT> 	 	Keep extending until number of extended contigs meet threshold, default 1

	--max_memory_TRINITY <INT> 		Maximum memory for TRINITY, default 20
	
	--threads <INT> 			Number of threads, default 8
   
	--threadsForSamSort <INT> 		Number of threads for samtools sort, default 2
	
	--memForSamSort <INT>			Memory for samtools sort, default 768M
   
	--complete_cleanup 			Complete cleanup of intermediate files, default partial cleanup
	
	--nochange_header 			Don't change header of final improved assembly file, default change header
	
	--cdHitEST <INT> 			0 for no CD-HIT-EST, 1 for CD-HIT-EST only in the start, 2 for CD-HIT-EST after every iteration, default 1
	
	--mapping_quality 			Minimum mapping quality to filter Bam file, default 20
	
	--sc_support_TH 			Softclips support for a position for Incomplete and fragmented contigs and false chimera process, default 0.75
	
	--min_SCs_reads <INT> 		  	Minimum number of softclips reads for extension of Incomplete contigs, default 3
	
	--sc_cons_len_TH <INT> 		 	Minimum length of consensus sequence generated by softclips to BLAST for fragmented contigs Identification, default 10
	
	--discard_contig_corner_len_TH <INT>  	Minimum number of bases allowed to discard from the corner of contigs while finding overlapped edges, default 10
	
	--max_allowed_gaps <INT> 		Maximum gaps allowed to find overlaps between the corners of the two contigs, default 0
	
	--edge_boundary <INT>			Terminal region for mapped reads to be considered for un-mapped and distantly mapped reads, default 2x of read length
	
	--min_unmapped_reads <INT> 		Minimum number of unmapped reads to generate CAP3 assembly for incomplete contigs extension, default 5
	
	--min_distMapped_reads <INT> 	 	Minimum number of distantly mapped reads to consider as read island for merging of fragmented contigs, default 3
	
	--contig_boundary <true> 		Check overlap between two fragmented contigs within 5% read length of the contig boundary, default read island boundary
	
	--min_allowed_unmapped_ext <INT> 	Minimum number of bases to be considered for valid extension using unmapped reads, default 50% of read length
	
	--one_side_allowed_SCs_RI <INT> 	Maximum % of softclips of the total read length allowed at one side of read to consider it for Read Island, default 25
	
	--each_side_allowed_SCs_RI <INT> 	Maximum % of softclips of the total read length allowed at both sides of read to consider it for Read Island, default 12
	
	--sc_start_pos_from_terminus <INT> 	For terminal softclips extraction, define starting position to consider, default 25
	
	--min_overlap_TH <INT> 		 	Minimum length of overlapped sequence to consider for merging fragmented contigs and CAP3 assemblies, default 20 bases
	
	--win_size <INT> 			Window size to detect gradual coverage change, default 2X read_length
	
	--win_diff_TH <INT> 			Maximum average coverage change threshold between two consecutive windows, default off
	
	--coverage_drop_TH <INT> 		Maximum coverage ratio between two consecutive positions to process for false chimera identification process, default 0.2
	
	--st_end_boundary <INT> 		Consider coverage change within start and end boundary of the contig, default one and half of read length 
	
	--blast_score_TH <INT> 		 	BLAST hit identity and coverage score for mis-assembly/false chimera, default 90
	
	--ignore_short_seq <INT> 		Minimum length of the contig to remove from the assembly, default 200
	
	--Insert_Ns <INT> 			Number of Ns to insert between contig and its CAP3 assembly in the absence of overlap, default 5
	
	Terminate ROAST process:		to stop ROAST properly before completion of default iterations place empty file named 'stop.txt' in the folder 'intermediate_Improved_assemblies'.

   
Output:
Improved fasta files from all iterations named as Name_SuperIterationNumber-MiniIterationNumber.fasta e.g. Name_1-1.fasta for first improvement in intermediate_improved_assemblies folder.

Final improved file named final_improved.fasta, summary file and removed sequences can be found in final_improved assembly folder.

Along with that, the log data for each improvement at each step is present in log folder with file name after each super and mini iterations.
Time log shows time taken by each step.

Run statistics display on screen.
