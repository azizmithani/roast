/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   byRI.h
 * Author: madiha
 *
 * Created on September 17, 2018, 1:53 PM
 */

#ifndef BYRI_H
#define BYRI_H
#include<map>
#include <string>
#include <set>
#include <vector>
#include <list>
#include <fstream> 
using namespace std;


class byRI {
public:
    byRI();
    byRI(const byRI& orig);
    virtual ~byRI();
    int read_island(string bam_in, string in_file,string outfile, string &log,  std::ofstream &log_time);
    string overlap_merge(string RI_seq, string MI_seq, int RI_start, int RI_end, int MI_start, int MI_end);
    //string overlap_merge(string str1, string str2);//, int str1_pos);
    int  merge_RIs(map<string, string> &final_merged_frags, map<string, vector <string> > &id_map, map<string, string> &merged_fragments, string contig_RI, string contig_MI, set<string> &dir_update, string RI_seq, string MI_seq, int RI_start, int RI_end, int MI_start, int MI_end, bool RI_dir, bool MI_dir, string header_fasta, string idd2, map<string, string> ::iterator &it_ri, ofstream &log_new, ofstream &mergedRI, set <string> &merged_check) ;
    void filter_merge(string merged_file, string merged_sam_file,string extended_file);
    void extend_unmappedBam(string bam_file, string in_file, string &out_file, string &log);
    

    
private:

};

#endif /* BYRI_H */

