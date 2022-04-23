/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
//#ifndef BYRI_H
#define BYRI_H
#include<map>
#include <string>
#include <vector>
#include <list>
#include <set>
#include "global.h"
using namespace std;

/*parameters about input data and filtering*/

int avg_IS = 0; // initialize

int read_length = 75;//100; // got value from average_insertsize() from utils which is being called in mini iteration

map<string, string> AllFasta_data;

int mapping_quality_TH = 20;

int base_quality = 20;

map<string, map<int, int> > coverage;

int contig_count = 0;

int max_coverage = 0;

vector<string> unmapped_contig;

set <string> extended_byCAP3;

map<string, list<int> > repeat_contigs;

/*Tags*/

string tool_name_flag = "tr";

string ROASTcap3_left_tag = "CAP3LT";

string ROASTcap3_right_tag = "CAP3RT";

string ROAST_LE_tag = "trLE";

string ROAST_RE_tag = "trRE";

string ROAST_BE_tag = "trBE";

string ROAST_update_tag = "trU";

string ROAST_newconitg1 = "a0";

string ROAST_newconitg2 = "b0";

string stop_LCAP3_ext = "CAP3LX";

string stop_RCAP3_ext = "CAP3RX";

/*SCs and RI parameters*/

float SC_support_TH = 75; //0.75 * 100;

//int min_sc_len = 5;      // min number of softclips from a read to be consider for further process

int SC_cons_len_TH = 10; // to consider consensus SC length for blast in fragment identification

int min_sc_reads = 3; //min number of soft clipped reads for extension # first it was 5, changed it for cap3 assembly 

float one_side_allowed_SCs_RI = 25; //% SCs percentage of read_length allowed in determining distantly mapped reads in finding broken by RI process 

float each_side_allowed_SCs_RI = 12; //% SCs percentage of read_length allowed in determining distantly mapped reads in finding broken by RI process

int min_unmapped_reads_CAP3 = 3;

int min_distMapped_reads = 5;

bool contig_boundary = false;

int min_CAP3_ext = 50; // 50% of read length

list<string>  cap3_left_list(contig_count), cap3_right_list(contig_count), split_contigs(contig_count);

multimap <string, int> merged_contigs;

int left_edge_boundary = read_length; //100;

int right_edge_boundary = left_edge_boundary * 2; //200;

int max_header_size = 1000; 

/*mis-assembly chimera parameters*/

int min_SCs_MCprocess = 25;

int max_SCs_oppositeSide_MCprocess = 3;

int sc_pos_from_corner = 25;

int consecutive_missAssembled_pos_dist = read_length;//  100;// 

int win_size = read_length * 2;//200; // bases

float coverage_drop = 0.2; //0.5

float avg_cov_th = 50;

float win_diff_TH = 0;// 0.2; // 0.25;
        
int st_end_th =  read_length;// 150; //1 n half of read length

int indel_TH = 10; // change TH to 100 for exempting exons

int mis_assembly_chimera_par_len = 75;

/*Parameters common for all*/

int blast_score_TH = 90; //0.80; // for mis-assembly/false chimera process

int discard_contig_cor_len = 10;

int min_overlap_TH  = 25; //20;

int max_allowed_gaps = 0;

/*Output parameters*/

int Ns = 5;

int min_extended_contigs = 1;

int improvment_TH = 1;   

int inner_itr_TH = 100;//30;

int outer_itr_TH = 100;//5;

int inner_iteration = 1, outer_iteration = 1,  ignore_short_seq = 200;

string delimiter = "\t";

int extended_incomplete_contigs = 0;

int merged_fragmented_contigs = 0;

int fixed_chimeric_contigs = 0;

int updatad_local_misAssemblies = 0;

string tool_name = "ROAST"; // Reference free Optimization of Assembled Supertranscriptomes

bool generate_assembly = false;

int cdhitest = 1; // 0:no, 1: only in the start, 2: at the end of every iteration

string out_directory = "";

bool complete_cleanup = false;

bool transRate_score = false;

bool change_header = true;

int temperory_contigs = 0, unmappded_count= 0 ;

string exe_path;

string path_inter;

/*Threads*/
string max_memory_TRINITY = "20";

int allowed_threads = 5;
string threads = "5";

int threadsForSamSort = 2;
string threadsForSamSortS = "2";

int MemForSamSort = 786;
string MemForSamSortS = "786";

int embed_count = 0;
