/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   extend_fasta.h
 * Author: madiha
 *
 * Created on November 5, 2019, 5:38 PM
 */

#ifndef BYSOFTCLIP_H
#define BYSOFTCLIP_H
#include<map>
#include <string>
#include <vector>
#include <list>
#include <fstream> 
#include <set>   
#include <map>
using namespace std;

class bySoftclip {
public:
    bySoftclip();
    bySoftclip(const bySoftclip& orig);
    virtual ~bySoftclip();
    bool extend_bySoftclip(string fastafile, string new_assembly, string bam_file, string iteration_log);
    int merge_fragments(map<string, string> &final_merged_frags, map<string, string> &merged_fragments, string contig_1, int contig1_length, int contig2_length, string contig_2, bool extended_side, string updated_seq, string contig2_seq, int contig1_st, int contig1_end, int contig2_st, int contig2_end, bool strand_contig2, string header_fasta, string idd2, ofstream &log, set<string> &merged_contigs2, map<string, bool> &ID_SCside, bool &fragments_merged);
    vector<string> extract_maxSC_block(vector<int> SC_base_pos, vector<string> SC_seq);
    int SC_start_pos(vector<int> SC_base_pos);
    string consensus_seq(vector<string>softclip_list, string flag); 
    int broken_by_blast(string new_assembly, string iteration_log, string extended_assembly);
    string overlapmerge_left(string str1, string str2, int str1_st, int str1_end, int str2_st);
    string overlapmerge_right(string str1, string str2, int str1_end, int str2_st, int str2_end) ;
    int fix_fragmentedBySCs(string fastafile, string new_assembly, string iteration_log); 
    //string overlapmerge_left(string str1, string str2, int str1_st, int str2_st);
    //string overlapmerge_right(string str1, string str2, int str1_end, int str2_end);
    char base_code(char base1, char base2);
    
      
private:

};

#endif /* EXTEND_FASTA_H */



