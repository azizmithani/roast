/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   global.h
 * Author: madiha
 *
 * Created on August 8, 2019, 2:05 PM
 */

#ifndef GLOBAL_H
#define GLOBAL_H
#include <map>
#include <set>
#include <string>
#include <vector>
#include <list>
using namespace std;

extern int avg_IS;

extern int read_length;

extern int mapping_quality_TH;

extern int base_quality;

extern map<string, map<int, int> > coverage;

extern vector<string> unmapped_contig;

extern set<string> extended_byCAP3;

extern int contig_count, max_coverage;

extern map<string, list<int> > repeat_contigs;

extern float SC_support_TH;

extern int SC_cons_len_TH;

extern int min_sc_len, unmappded_count;

extern list<string> cap3_left_list, cap3_right_list , split_contigs;

extern multimap <string, int> merged_contigs;

extern string path_inter;

extern string delimiter;

extern int improvment_TH;

extern int min_unmapped_reads_CAP3;

extern int min_distMapped_reads;

extern bool contig_boundary;

extern float one_side_allowed_SCs_RI;

extern float each_side_allowed_SCs_RI; 

extern string tool_name_flag;

extern string ROASTcap3_left_tag;

extern string ROASTcap3_right_tag;

extern string ROAST_LE_tag;

extern string ROAST_RE_tag;

extern string ROAST_BE_tag;

extern string ROAST_update_tag;

extern string ROAST_newconitg1;

extern string ROAST_newconitg2;

extern string stop_LCAP3_ext;

extern string stop_RCAP3_ext;

extern int min_CAP3_ext;

extern int min_SCs_MCprocess;

extern int max_SCs_oppositeSide_MCprocess;

extern int extended_incomplete_contigs;

extern int merged_fragmented_contigs;

extern int fixed_chimeric_contigs;

extern int updatad_local_misAssemblies;

extern int min_sc_reads;

extern int sc_pos_from_corner;

extern int left_edge_boundary;

extern int right_edge_boundary;

extern int max_header_size;

extern int discard_contig_cor_len;

extern int min_overlap_TH;

extern int consecutive_missAssembled_pos_dist;// 

extern int blast_score_TH; 

extern int win_size;

extern float coverage_drop, avg_cov_th, win_diff_TH; // try3 70;

extern int st_end_th; //1 n half of read length

extern int indel_TH;

extern int mis_assembly_chimera_par_len;

extern int Ns;

extern int min_extended_contigs;

extern int min_overlap_TH ;

extern int max_allowed_gaps;

extern string tool_name;

extern bool generate_assembly;

extern string TransRate_Home;

extern int inner_itr_TH ;

extern int outer_itr_TH ;

//extern bool basic_cleanup ;
extern bool complete_cleanup;

extern bool transRate_score ;

extern bool change_header ;

extern int cdhitest;

extern int inner_iteration;

extern int outer_iteration;

extern int temperory_contigs;

extern int ignore_short_seq;

extern string max_memory_TRINITY;

extern int allowed_threads;
extern string threads;

extern int threadsForSamSort;
extern string threadsForSamSortS;

extern int MemForSamSort;
extern string MemForSamSortS;

extern string out_directory;

extern string exe_path;

extern map<string, string> AllFasta_data;

extern int embed_count;

#endif /* GLOBAL_H */
