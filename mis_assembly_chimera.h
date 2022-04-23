/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   mis_assembly_chimera.h
 * Author: madiha
 *
 * Created on December 29, 2019, 10:11 PM
 */

#ifndef MIS_ASSEMBLY_CHIMERA_H
#define MIS_ASSEMBLY_CHIMERA_H
#include <map>
#include <string>
#include <vector>
#include <list>
using namespace std;

class mis_assembly_chimera {
public:
    mis_assembly_chimera();
    mis_assembly_chimera(const mis_assembly_chimera& orig);
    virtual ~mis_assembly_chimera();
  

   int extract_chimera_mis_assembly(string bam_in, string coverage_file_hisat, string coverage_file_bowtie, string fastafile, string processed_fastafile, string cufflink_file, string mis_assembly_chimera_log);
    
   string chimera_type(int pos, string contig_seq);
   
   vector<int> chimera_cov1(string contig_name, vector<int> coverage, int cov_th);
   
   vector< int>gradual_change(string contig_name, vector<int> coverage, int cov_th);
   
   vector<int> single_SC_pattern(string contig_name, vector<int>left_SC_base_pos, vector<string>left_SC_seq, vector<int>right_SC_base_pos, vector<string>right_SC_seq, vector<int> cov, int cov_th);
   
   vector<int> compare_singleSC(vector<int> singleSC_list, list<int> ex_bo, list<int> sc_pat);
   
   vector<int> compare_gradual_cov_change(vector<int>gradual_dip_chimera, list<int> ex_bo, vector<int> sc_pat); 
  
   vector<int> sc_pattern(string contig_name,string contig_seq, multimap<int, string> left_SC_pos_seq, multimap<int, string> right_SC_pos_seq, vector<int>cov, string &output_file, string PATH, string bam_in);
   
   vector<int> filter_singleSC(vector<int> single_SC, map <int, string> unmappedMate_reads, map<int, string> distMapped_reads);
   
   bool filter_missingSeq(string contig_name, int last_rightSC, string bam_file);
   
   string perform_blast(string query_seq, string subject1_seq, string subject2_seq);
  
   string perform_blast(string query_seq, string subject1_seq);
 
   string overlapSC(string RSC, string LSC);
 
   int filterPos_processFasta(string fastafile , string processed_fastafile, string identified_posFile, string mis_assembly_chimera_log, map<string, list<int> > &exon_boundary);
  
   void get_problemData(vector <string> same_contig, vector<int> &positions,vector<int> &problem_type);
  
   void update_positions(vector<int> &all_positions, int start, int end, int update_by, bool to_push);
  
   bool compare_positions(vector<int> all_positions, vector<int> problems, int p_start, int p_end, int problem);
 
   bool compare_positions(vector<int> all_positions, vector<int> problems, int p_start, int p_end, int insert_at, int problem );
  
   bool compare_positions(vector<int> all_positions, vector<int> problems, int pos , int problem);
private:
    


};

#endif /* MIS_ASSEMBLY_CHIMERA_H */
