/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   extend_fasta.cpp
 * Author: madiha
 * 
 * Created on November 5, 2019, 3:38 PM
 */

#include "bySoftclip.h"
#include "ROAST_extendContigs_SCs.h"
#include "ROAST_mergeContigs_SCs.h"
//#include "ROAST_mergeContigs_SCs_R2.h"
#include "utils.h"
#include "global.h"
#include "alignment.h"
#include <sstream>
#include <iostream>
#include <fstream>                                                                                                                                                                                         
#include <algorithm>
#include <vector>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <cctype>    
#include <climits>
#include <set>
#include <math.h>   
#include <map>
#include <list>
#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <boost/regex.hpp>
#include <boost/filesystem/operations.hpp>

using namespace BamTools;
using namespace std;
using namespace boost;

bySoftclip::bySoftclip() {
}

bySoftclip::bySoftclip(const bySoftclip& orig) {
}

bySoftclip::~bySoftclip() {
}

struct fragments_data {
    string contig1;
    int contig1_start;
    int contig1_end;
    bool contig1_strand; // reverse = true, forward = false //string strand_mate1;
    bool isLSC;

    int BLAST_score;
    string SC_seq;

    string contig2;
    int contig2_length;
    int contig2_start;
    int contig2_end;
    bool contig2_strand; // reverse = true, forward = false //string strand_mate1;
    bool flag;


};

struct merged_info {
    string updated_id;
    string updated_seq;
    string updated_strand; // forward or RC for merging

};
utils utils;

bool bySoftclip::extend_bySoftclip(string fastafile, string new_assembly, string bam_file, string iteration_log) {

    /*do alignment first and calculate coverage*/
   int contig_counts = 1, count = 1, thread_count = 1;
    typedef multimap<string, string > mmap;
    mmap SC_id;
    typedef pair<string, string> left_right;

    int extended_count = 0;
    bool extend;

    vector<string> temp_vec;
    std::stringstream ct;
    ct << thread_count;

    ifstream fasta_file;
    fasta_file.open(fastafile.c_str());

    string ID_file = path_inter + "/IDs_" + ct.str() + ".bed";
    string subfasta_file = path_inter + "/thread_" + ct.str() + ".fasta"; // write for each 

   // ct.str(string());

    ofstream IDfile;
    IDfile.open(ID_file.c_str());

    ofstream sub_fastaFile;
    sub_fastaFile.open(subfasta_file.c_str());

    //  cout << contig_count << endl;
    //contig_count = 48295;
    int sub_fasta;
    if (contig_count > allowed_threads && allowed_threads > 1)
        sub_fasta = round(contig_count / allowed_threads);
    else
        sub_fasta = contig_count;

    string entry;
    string seq, ID;

    string ROAST_extendContigs_SCs;
    stringstream P1, P2, P3, P4;
    P1 << SC_support_TH;
    P2 << min_sc_reads;
    P3 << sc_pos_from_corner;
    P4 << exe_path;
    //write contigs names in n number of files based on number of threads
    map<string, string>::iterator it_fasta;

    while (!(getline(fasta_file, entry).eof())) {

        if (entry[0] == '>') { // get contig IDs for currently generated fasta file to extract respective bam region
            ID = entry.substr(1, entry.size());

            if (thread_count <= allowed_threads - 1) {

                if (contig_counts <= sub_fasta) { // write n numbers of contigs in each fasta file
                    it_fasta = AllFasta_data.find(ID);
                    seq = it_fasta->second;
                    sub_fastaFile << entry << endl << seq << endl;
                    IDfile << ID << "\t" << "0" << "\t" << seq.size() << endl; // ID start_pos end_pos
                    contig_counts++;
                    //stringstream seq;
                } else { // close current file and open new file for n contigs again
                    sub_fastaFile.close();
                    IDfile.close();

                    ROAST_extendContigs_SCs = exe_path + "ROAST_extendContigs_SCs " + subfasta_file + " " + ID_file + " " + bam_file + " " + ct.str() + " " + P1.str() + " " + tool_name_flag + " " + P2.str() + " " + P3.str() + " " + P4.str() + " &";
                    std::system(ROAST_extendContigs_SCs.c_str());
        
                    thread_count++;
                    ct.str(string()); // clear for previous thread
                    ct << thread_count;

                    //sub_bam = path_inter + "/thread_" + ct.str() + ".bam"; // write for each fasta file 
                    ID_file = path_inter + "/IDs_" + ct.str() + ".bed";
                    subfasta_file = path_inter + "/thread_" + ct.str() + ".fasta";

                    sub_fastaFile.open(subfasta_file.c_str());
                    IDfile.open(ID_file.c_str());

                    it_fasta = AllFasta_data.find(ID);
                    seq = it_fasta->second;
                    sub_fastaFile << entry << endl << seq << endl;
                    IDfile << ID << "\t" << "0" << "\t" << seq.size() << endl; // ID start_pos end_pos

                    contig_counts = 2; // because first entry has already been added
                   // ct.str(string()); //clear();
                    
                }
            } else { // add remaining in last thread
                it_fasta = AllFasta_data.find(ID);
                seq = it_fasta->second;
                sub_fastaFile << entry << endl << seq << endl;
                IDfile << ID << "\t" << "0" << "\t" << seq.size() << endl; // ID start_pos end_pos

                contig_counts++;
                //stringstream seq;
            }
        }

    }

    sub_fastaFile.close();
    IDfile.close();
    // for last thread 15nov21
    ROAST_extendContigs_SCs = exe_path + "ROAST_extendContigs_SCs " + subfasta_file + " " + ID_file + " " + bam_file + " " + ct.str() + " " + P1.str() + " " + tool_name_flag + " " + P2.str() + " " + P3.str() + " " + P4.str() +" &";
    std::system(ROAST_extendContigs_SCs.c_str());

    // ct.str(string()); //clear();
    //call for all threads

    /*for (int i = 1; i <= allowed_threads; i++) {
        ct << i;
        //extend contigs for each bam file using Soft clips at terminal
        subfasta_file = path_inter + "/thread_" + ct.str() + ".fasta";
        ID_file = path_inter + "/IDs_" + ct.str() + ".bed";

        //ROAST_extendContigs_SCs tt;
        //tt.extend_bySoftclips(subfasta_file, ID_file, bam_file, ct.str(), SC_support_TH, tool_name_flag, min_sc_reads, sc_pos_from_corner);

        ROAST_extendContigs_SCs = exe_path + "ROAST_extendContigs_SCs " + subfasta_file + " " + ID_file + " " + bam_file + " " + " " + ct.str() + " " + P1.str() + " " + tool_name_flag + " " + P2.str() + " " + P3.str() + " &";
        std::system(ROAST_extendContigs_SCs.c_str());

        ct.str(string());
    }*/

    ct.str(string());
    P1.str(string());
    P2.str(string());
    P3.str(string());
    P4.str(string());

    int running_threads_count = allowed_threads;

    while (running_threads_count > 0) {
        for (int i = 1; i <= allowed_threads; i++) {

            ct << i;
            string thread_done = path_inter + "/done_" + ct.str() + ".txt";
            ct.str(string());

            if (boost::filesystem::exists(thread_done)) // does filePath actually exist?
            {
                //continue;             
                running_threads_count--;

            } else { // unless exist status for all files obtained
                // i = 1; // start from 1st file again to check status
                sleep(10); // 
                running_threads_count = allowed_threads;
                break;
            }

        }
    }

    std::size_t found = new_assembly.find_last_of("/\\");
    string extendfasta = new_assembly ;//new_assembly.substr(0, found) + "/ext.fasta"; // no merginf now so extend fasta is final improved file

    ofstream extended_fasta;
    extended_fasta.open(extendfasta.c_str());

    string merge_assemblies = "cat ";
    string merge_logs = "cat ";

    for (int i = 1; i <= allowed_threads; i++) {

        ct << i;
        string individual_extended_files = path_inter + "/ext_" + ct.str() + ".fasta";
        string individual_log_files = path_inter + "/iteration_log_" + ct.str() + ".txt";

        merge_assemblies = merge_assemblies + individual_extended_files + " ";
        merge_logs = merge_logs + individual_log_files + " ";
        ct.str(string());
    }
    ct.str(string());
    merge_assemblies = merge_assemblies + " > " + extendfasta;
    merge_logs = merge_logs + " > " + iteration_log;

    std::system(merge_assemblies.c_str());
    std::system(merge_logs.c_str());

    //remove all intermediate file of this process
    for (int i = 1; i <= allowed_threads; i++) {

        ct << i;
        string individual_extended_files = path_inter + "/ext_" + ct.str() + ".fasta";
        string individual_log_files = path_inter + "/iteration_log_" + ct.str() + ".txt";
        string sub_bam = path_inter + "/thread_" + ct.str() + ".bam";
        string subfasta_file = path_inter + "/thread_" + ct.str() + ".fasta";
        string sub_cov = path_inter + "/thread_" + ct.str() + ".cov";
        string ID_file = path_inter + "/IDs_" + ct.str() + ".bed";

        string extended_bySCs_file = path_inter + "/SC_id_" + ct.str() + ".txt";
        remove(extended_bySCs_file.c_str());

        string thread_done = path_inter + "/done_" + ct.str() + ".txt";
        ifstream done;
        done.open(thread_done.c_str());
        getline(done, entry);
        extended_count = extended_count + atoi(entry.c_str());
        done.close();

        remove(thread_done.c_str());
        remove(individual_extended_files.c_str());
        remove(individual_log_files.c_str());
        remove(sub_bam.c_str());
        remove(sub_cov.c_str());
        remove(subfasta_file.c_str());
        remove(ID_file.c_str());
        ct.str(string());

    }
   // ct.str(string());

    cout << "Extended incomplete contigs using soft-clipping at iteration " << outer_iteration << "-" << inner_iteration << ": " << extended_count << endl;
    extended_fasta.close();
    // utils.remove_file(cov_file);
    if (extended_count >= min_extended_contigs)
        extend = true;
    else
        extend = false;

    extended_incomplete_contigs = extended_incomplete_contigs + extended_count;


    /*Merge Frgamented contigs using SCs in threads*/
  /* string entry_ID;
      vector<string> temp;
      string command_ROAST_mergeContigs_SCs;

      for (int i = 1; i <= allowed_threads; i++) {
          ct << i;
          P1 << SC_cons_len_TH;
          P2 << min_overlap_TH;
          P3 << blast_score_TH;
          P4 << max_allowed_gaps;


          string extended_bySCs_file = path_inter + "/SC_id_" + ct.str() + ".txt";
          command_ROAST_mergeContigs_SCs = exe_path + "ROAST_mergeContigs_SCs " + ct.str() + " " + extendfasta + " " + path_inter + " " + P1.str() + " " + exe_path + " " + P2.str() + " " + P3.str() + " " + P4.str()  + " &";
          std::system(command_ROAST_mergeContigs_SCs.c_str());

          //tms.FindFragmentsBySCs( ct.str(), extendfasta, path_inter);
          ct.str(string());
          P1.str(string());
          P2.str(string());
          P3.str(string());
          P4.str(string());
      }

      running_threads_count = allowed_threads;

      while (running_threads_count > 0) {
          for (int i = 1; i <= allowed_threads; i++) {

              ct << i;
              string thread_done = path_inter + "/done" + ct.str() + ".txt";
              ct.str(string());

              if (boost::filesystem::exists(thread_done)) // does filePath actually exist?
              {
                  //continue;             
                  running_threads_count--;

              } else { // unless exist status for all files obtained
                  // i = 1; // start from 1st file again to check status
                  sleep(10); // 
                  running_threads_count = allowed_threads;
                  break;
              }

          }
      } 
    for (int i = 1; i <= allowed_threads; i++) {
        ct << i;
        string thread_done = path_inter + "/done" + ct.str() + ".txt";
        string extended_bySCs_file = path_inter + "/SC_id_" + ct.str() + ".txt";

        utils.remove_file(thread_done.c_str());
        utils.remove_file(extended_bySCs_file.c_str());
        ct.str(string());

    }
     broken_by_blast(new_assembly, iteration_log, extendfasta, merged_contigs);
    utils.remove_file(extendfasta); */

    return extend;

}

struct length {

    bool operator()(const string& a, const string& b) {
        return a.size() < b.size();
    }
};

int bySoftclip::broken_by_blast(string new_assembly, string iteration_log, string extended_assembly) {


    alignment alignment;
    AllFasta_data = alignment.extract_fasta_data(extended_assembly); // load whole extended assembly in global variable 

    ofstream mergedAssembly;
    mergedAssembly.open(new_assembly.c_str());
    ofstream log;

    log.open(iteration_log.c_str());//), fstream::app);
    vector<string> temp2, same_contig;
    set<string> merged_contigs2; // 2 because merged_contigs is for over all merged sequences, but to track merged contigs of this iteration we need separate  
    same_contig.reserve(10);
    //merged_contigs2.reserve(contig_count);

    typedef multimap<string, string > mmap;
    //ofstream broken_contig;

    typedef multimap<int, fragments_data> fragments;
    fragments broken_bySCs;

    vector <string> temp;

    string entry2, contig_name, contig_left_merged, contig1_seq, contig2_seq, contig_1, contig_2, SC_seq, idd1, idd2, contig2_seq_RC, ID, entry_ID, merge;
    int contig1_st, contig2_st, contig1_end, contig2_end, contig1_length, contig2_length, contig2_st_rc, contig2_end_rc;
    int count = 0, merged_count = 0;
    //float blast_score_TH = 80;

    fragments_data fragments_data;

    stringstream ct;
    for (int i = 1; i <= allowed_threads; i++) {

        ct << i;
        string filterd_frags = path_inter + "/filtered_frags_" + ct.str() + ".txt";
        ifstream filt;
        filt.open(filterd_frags.c_str());

        while (!(getline(filt, entry_ID).eof())) {

            temp.clear();
            utils.str_split(entry_ID, temp, delimiter);

            if (temp[0] == "left")
                fragments_data.isLSC = true;
            else
                fragments_data.isLSC = false; // right

            if (temp[4] == "+")
                fragments_data.contig2_strand = false; // +
            else
                fragments_data.contig2_strand = true; // -

            fragments_data.contig1_strand = false; // as SCs are extended on positive strands already
            
            fragments_data.contig1 = temp[1];
            fragments_data.contig1_start = atoi(temp[5].c_str());
            fragments_data.contig1_end = atoi(temp[6].c_str());

            fragments_data.contig2 = temp[2];
            fragments_data.BLAST_score = atoi(temp[3].c_str());

            fragments_data.contig2_start = atoi(temp[8].c_str());
            fragments_data.contig2_end = atoi(temp[9].c_str());
            fragments_data.contig2_length = atoi(temp[7].c_str());
            fragments_data.SC_seq = temp[11];

            fragments_data.flag = true;

            broken_bySCs.insert(make_pair(fragments_data.BLAST_score, fragments_data));

        }

        filt.close();
        remove(filterd_frags.c_str());
        ct.str(string());
    }

    // sorted map based on max BLAST score per following rules
    // saved in ascending sorted order by read count
    // for same outer key it get the first one, otherwise it sort based on key of inner map
    // and if inner keys are same then it sort based on outer keys..

    typedef map<string, string> merged_fragment; // entries of this vectors are contig1 + header_fasta, contig2 + header_fasta
    merged_fragment merged_fragments;

    typedef map<string, bool > id_SCside; // header_fasta + contig1 + contig2
    id_SCside ID_SCside;

    typedef map<string, string> final_merged_frag; // header_fasta + new merged sequence
    final_merged_frag final_merged_frags;
    set<string> dir_update;

    fragments::reverse_iterator it_SCs;
    fragments::iterator it_data;
    string new_seq;
    map<string, string>::iterator it_fasta;
    bool fragments_merged;
    //cout << "loop over whole map, size of broken_bySCs: " << broken_bySCs.size() << endl;


    for (it_SCs = broken_bySCs.rbegin(); it_SCs != broken_bySCs.rend(); it_SCs++) { //map<sting, map <int, RI_data> >
        // cout << it_SCs->first << " " << it_SCs->second.contig1 << endl;
        string header_fasta;
        stringstream ss;

        contig_1 = it_SCs->second.contig1;
        contig_2 = it_SCs->second.contig2;

       // cout << "contig 1: " << contig_1 << " contig 2: " << contig_2 << endl;

        bool strand_contig1 = it_SCs->second.contig1_strand; // RSC true LSC true 
        bool strand_contig2 = it_SCs->second.contig2_strand; // - true + false
        bool isLSC = it_SCs->second.isLSC;

        contig2_st = it_SCs->second.contig2_start;
        contig2_end = it_SCs->second.contig2_end;

        contig2_length = it_SCs->second.contig2_length;
        SC_seq = it_SCs->second.SC_seq;
        bool flag = it_SCs->second.flag;

        it_fasta = AllFasta_data.find(contig_1);
        contig1_seq = it_fasta->second;

        it_fasta = AllFasta_data.find(contig_2);
        contig2_seq = it_fasta->second;

        contig1_st = it_SCs->second.contig1_start;
        contig1_end = it_SCs->second.contig1_end;
            
        //std::size_t contig1_st = contig1_seq.find(SC_seq);
       // contig1_end = contig1_seq.size() - 1;
        contig1_length = contig1_seq.size();

        /*if(!isLSC){ //if (it_SCs->second.contig1_start != 0) { // right // map because start and end positions of query are for SC sequences only
            contig1_st = contig1_seq.size() - it_SCs->second.contig1_end;
            contig1_end = contig1_seq.size() - it_SCs->second.contig1_start + 1;
        }*/
        // else for Left SCs same start and end position
        idd1 = contig_1;
        idd2 = contig_2; //.substr(10, contig_2.size());  

        size_t update_ind = contig_1.find(tool_name_flag);

        if (update_ind != std::string::npos) { // [] found remove it and update 
            idd1 = contig_1.substr(0, update_ind - 1);
        }

        update_ind = contig_2.find(tool_name_flag);

        if (update_ind != std::string::npos) { // [] found remove it
            idd2 = contig_2.substr(0, update_ind - 1);
        }

        // cout << contig_1 << " and " << contig_2 << endl;
        if (flag) {
            string ch1 = contig_1;
            string ch2 = contig_2;

            string ch3 = contig_2;
            string ch4 = contig_1;

            size_t chimera11 = contig_1.rfind(ROAST_newconitg1); // because multi-gene Chimera has _a and _b at end of t$
            size_t chimera12 = contig_2.rfind(ROAST_newconitg2);

            size_t chimera21 = contig_2.rfind(ROAST_newconitg1); // because multi-gene Chimera has _a and _b at end of t$
            size_t chimera22 = contig_1.rfind(ROAST_newconitg2);

            //   if (contig_RI != contig_MI) { have already check in the start
            if (chimera11 != string::npos)
                ch1 = contig_1.substr(0, chimera11 - 1);

            if (chimera12 != string::npos)
                ch2 = contig_2.substr(0, chimera12 - 1);

            if (chimera21 != string::npos)
                ch3 = contig_2.substr(0, chimera21 - 1);

            if (chimera22 != string::npos)
                ch4 = contig_1.substr(0, chimera22 - 1);


            if ((ch1 != ch2) && (ch4 != ch3)) { //either both contigs are not same or if they are chimeras split by _a, _b then their original contig ids should not same before _a and _b
                // merge them
                // find read contig if read strand is same change flag to false
                // find mate contig (which appears as read now)if in same direction as just dealt here make flag false
                // we also need read strand information and only change read contig direction to make it trackable


                //check if contig_RI or contig_MI is already updated? use their updated name, seq, and strands

                map<string, string> ::iterator it_sc1;
                it_sc1 = merged_fragments.find(contig_1);

                map<string, string> ::iterator it_sc2;
                it_sc2 = merged_fragments.find(contig_2);

                if (it_sc1 != merged_fragments.end()) {
                    //      cout << "Read contig is repeating " << endl;
                    header_fasta = it_sc1->second; //
                    final_merged_frag::iterator it_update_seq;
                    id_SCside::iterator it_side;

                    it_update_seq = final_merged_frags.find(header_fasta);
                    string updated_seq = it_update_seq->second;

                    it_side = ID_SCside.find(header_fasta); // LSC or RSC?
                    bool extended_side = it_side->second;
                    
                    merged_count = merged_count + merge_fragments(final_merged_frags, merged_fragments, contig_1, contig1_length, contig2_length, contig_2, extended_side, updated_seq, contig2_seq, contig1_st, contig1_end, contig2_st, contig2_end, strand_contig2, header_fasta, idd2, log, merged_contigs2, ID_SCside, fragments_merged);

                } else if (it_sc2 != merged_fragments.end()) {

                    continue;

                } else { // both not already written
                    //          cout << "new read and mate contigs" << endl;
                    if (isLSC) { // LSC 

                        if (!strand_contig2) { // +
                            // cout << "LSC +" << endl;
                      //      cout << idd2 << " " << idd1 << endl;
                            merge = overlapmerge_left(contig2_seq, contig1_seq, contig2_st, contig2_end, contig1_st); //0);
                            //merge_left = merge;

                            if (merge != "") {

                                ss << idd2 << "_and_" << idd1; // add [m] at the end or in new_header
                                merged_count++;
                                header_fasta = ss.str();
                                
                                if(header_fasta.length() > max_header_size) //27dec 21  to avoid error in cufflink due to long ID size
                                    header_fasta = idd2;             
                                
                                //28oct21 overlap boundry ectraction to avoid split of same position in multigene chimera
                                int start_overlap = contig2_st; // start of overlap
                                int end_overlap = contig2_end; // end of overlap
                                merged_contigs.insert(make_pair(header_fasta, start_overlap));
                                merged_contigs.insert(make_pair(header_fasta, end_overlap));

                                merged_contigs2.insert(contig_2);
                                ss.str(string());
                                                              
                                log << contig_1 << "\t" << contig_2 << "\t" << "5' side" << "\t " << "+strand" << "\t" << start_overlap << "-" << end_overlap << endl;
                            }
                        } else { // if (strand_contig2) { // - / take reverse complement
                            //  cout << " LSC -" << endl;
                            // dir_update.insert(strand_contig2);
                            contig2_seq_RC = utils.Rcomplement(contig2_seq);
                            contig2_st_rc = contig2_length - contig2_end;
                            contig2_end_rc = contig2_length - contig2_st;

                        //    cout << idd2 << " " << idd1 << endl;
                            merge = overlapmerge_left(contig2_seq_RC, contig1_seq, contig2_st_rc, contig2_end_rc,  contig1_st); // contig_st = 0
                            //merge_left = merge;

                            if (merge != "") {

                                ss << idd2 << "_and_" << idd1; // add [m] at the end or in new_header
                                merged_count++;
                                header_fasta = ss.str();
                                ss.str(string());
                                if (header_fasta.length() > max_header_size) //27dec 21  to avoid error in cufflink due to long ID size
                                    header_fasta = idd2;
                                
                                //28oct21 overlap boundry ectraction to avoid split of same position in multigene chimera
                                int start_overlap = contig2_st_rc; // start of overlap
                                int end_overlap = contig2_end_rc; // end of overlap
                                merged_contigs.insert(make_pair(header_fasta, start_overlap));
                                merged_contigs.insert(make_pair(header_fasta, end_overlap));

                                // merged_contigs2.push_back(contig_1);
                                merged_contigs2.insert(contig_2);

                                log << contig_1 << "\t" << contig_2 << "\t" << "5' side" << "\t " << "-strand" << "\t" << start_overlap << "-" << end_overlap << endl;

                            }

                        }
                    } else { // if (!LSC) { //RSC
                        // cout << "RSC +" << endl;
                        if (!strand_contig2) { // + 

                            //cout << idd1 << " " << idd2 << endl;
                            merge = overlapmerge_right(contig1_seq, contig2_seq, contig1_end, contig2_st, contig2_end);

                            if (merge != "") {

                                ss << idd1 << "_and_" << idd2; // << "_" << tool_name_flag << "_M_" << tool_name_flag; // "_merged"; // repeating things here, fixed*** // 
                                merged_count++;
                                header_fasta = ss.str();

                                if (header_fasta.length() > max_header_size) //27dec 21  to avoid error in cufflink due to long ID size
                                    header_fasta = idd1;

                                //28oct21 overlap boundry ectraction to avoid split of same position in multigene chimera
                                int start_overlap = contig1_st; // start of overlap
                                int end_overlap = contig1_end; // end of overlap
                                merged_contigs.insert(make_pair(header_fasta, start_overlap));
                                merged_contigs.insert(make_pair(header_fasta, end_overlap));
                                
                                //mergedAssembly << ">" << new_header << endl << merge << endl;
                                log << contig_1 << "\t" << contig_2 << "\t" << "3' side" << "\t " << "+strand" << "\t" << start_overlap << "-" << end_overlap << endl;

                                ss.str(string());
                            }

                        } else { // if (strand_contig2) { //-
                            // cout << "RSC - " << endl;
                            contig2_seq_RC = utils.Rcomplement(contig2_seq);
                            contig2_st_rc = contig2_length - contig2_end;
                            contig2_end_rc = contig2_length - contig2_st;
                           
                            //cout << idd1 << " " << idd2 << endl;
                            merge = overlapmerge_right(contig1_seq, contig2_seq_RC,  contig1_end, contig2_st_rc, contig2_end_rc);

                            if (merge != "") {

                                ss << idd1 << "_and_" << idd2; // << "_" << tool_name_flag << "_M_" << tool_name_flag;
                                merged_count++;
                                header_fasta = ss.str();

                                if (header_fasta.length() > max_header_size) //27dec 21  to avoid error in cufflink due to long ID size
                                    header_fasta = idd1;

                                //28oct21 overlap boundry ectraction to avoid split of same position in multigene chimera
                                int start_overlap = contig1_st; // start of overlap
                                int end_overlap = contig1_end; // end of overlap
                                merged_contigs.insert(make_pair(header_fasta, start_overlap));
                                merged_contigs.insert(make_pair(header_fasta, end_overlap));
                                
                                ss.str(string());
                                merged_contigs2.insert(contig_1);
                                merged_contigs2.insert(contig_2);
                                
                                                                // mergedAssembly << ">" << new_header << endl << merge << endl;
                                log << contig_1 << "\t" << contig_2 << "\t" << "3' side" << "\t " << "-strand" << "\t" << start_overlap << "-" << end_overlap << endl;
                            }
                        }
                    }
                    if (!(merge.empty()) && (header_fasta != "")) { //fill all four data structures
                        //  cout << "Fill out all DS" << endl;
                        merged_fragments.insert(make_pair(contig_1, header_fasta));
                        merged_fragments.insert(make_pair(contig_2, header_fasta));

                        vector <string> temp; // add new and vector of old ids in id map
                        temp.push_back(contig_1);
                        temp.push_back(contig_2);
                        ID_SCside.insert(make_pair(header_fasta, isLSC));
                        temp.clear();

                        final_merged_frags.insert(make_pair(header_fasta, merge)); // add new id and new sequence

                        merged_contigs2.insert(contig_1);
                        merged_contigs2.insert(contig_2);
                        fragments_merged = true;
                    }
                }
                // for current and all new entries if same contig RI appear with same strand? ignore it  
                if (fragments_merged) {

                    for (it_data = broken_bySCs.begin(); it_data != broken_bySCs.end(); it_data++) { // treat anchor contig only (both in RSC LSC)

                        if (it_data->second.contig1 == contig_2 || it_data->second.contig2 == contig_2)
                            it_data->second.flag = false;
                    }
                }
                fragments_merged = false;
                header_fasta = "";
            } //chimera check
        } // flag

        // find make all flags false here containing same contig_RI and direction and contig_MI and direction in broken_byRI map<sting, map <int, RI_data> >
    }

    for (final_merged_frag::iterator it_merged = final_merged_frags.begin(); it_merged != final_merged_frags.end(); it_merged++) {

        mergedAssembly << ">" << it_merged->first << endl << it_merged->second << endl;
        //  cout << ">" << it_merged->first << endl << it_merged->second << endl;
    }


    string ID_temp, entry_fasta;
    ifstream fastafile;
    fastafile.open(extended_assembly.c_str());

    while (!(getline(fastafile, entry_fasta).eof())) {

        if (entry_fasta[0] == '>') {
            //temp = utils.split(entry_fasta, " ");
            //temp.clear();
            utils.str_split(entry_fasta, temp, " ");
            ID_temp = entry_fasta;
            ID = entry_fasta.substr(1, entry_fasta.size());

        } else {

            if (find(merged_contigs2.begin(), merged_contigs2.end(), ID) == merged_contigs2.end()) { //current contig is not being merged

                mergedAssembly << ID_temp << endl << entry_fasta << endl;

            } else { // when ID found in list of merged sequences
                continue;
            }
        }
    }

    log.close();
    mergedAssembly.close();

    cout << "Merged incomplete contigs using Soft-clipping at iteration " << outer_iteration << "-" << inner_iteration << ": " << merged_count << endl;
    merged_fragmented_contigs = merged_fragmented_contigs + merged_count;
   // exit(0);
    return merged_count;
}


int bySoftclip::fix_fragmentedBySCs(string fastafile, string new_assembly, string iteration_log) {

    int contig_counts = 1, count = 1, thread_count = 1;
    std::stringstream ct;
    ct << thread_count;

    string ID_file = path_inter + "/IDs_" + ct.str() + ".txt";


    ofstream IDfile;
    IDfile.open(ID_file.c_str());
    
    int running_threads_count = allowed_threads;
        /*Merge Frgamented contigs using SCs in threads*/
    string entry_ID, ID;
    vector<string> temp;
    set <string> IDs;
    string command_ROAST_mergeContigs_SCs;

    ifstream in_file;
    in_file.open(fastafile.c_str());


     
    while (!getline(in_file, entry_ID).eof()) { //save all IDs of fastafile
        if (entry_ID[0] == '>') {
            ID = entry_ID.substr(1, entry_ID.size() - 1);

            if (ID.find("_a0") == string::npos && ID.find("_b0") == string::npos) {
                if (ID.find(ROAST_LE_tag) != string::npos || ID.find(ROAST_BE_tag) != string::npos || ID.find(ROAST_RE_tag) != string::npos) {
                    //  cout << ID <<endl;
                    IDs.insert(ID);
                }

            }
        }
    }
    in_file.close();
    
    int sub_IDs;
    
    if (IDs.size() > allowed_threads && allowed_threads > 1)
        sub_IDs = round(IDs.size() / allowed_threads);
    else
        sub_IDs = IDs.size();

    stringstream P1, P2, P3, P4;
    P1 << SC_cons_len_TH;
    P2 << min_overlap_TH;
    P3 << blast_score_TH;
    P4 << max_allowed_gaps;
    
    for (set<string>::iterator it_ID = IDs.begin(); it_ID != IDs.end(); it_ID++) {
        if (thread_count <= allowed_threads - 1) {

            if (contig_counts <= sub_IDs) { // write n numbers of contigs in each fasta file

                IDfile << *it_ID << endl; // IDs
                contig_counts++;
                //stringstream seq;
            } else { // close current file and open new file for n contigs again
                IDfile.close();

                command_ROAST_mergeContigs_SCs = exe_path + "ROAST_mergeContigs_SCs " + ct.str() + " " + fastafile + " " + path_inter + " " + exe_path + " " + P1.str() + " " + P2.str() + " " + P3.str() + " " + P4.str() + " &";
                std::system(command_ROAST_mergeContigs_SCs.c_str());
                
                ct.str(string());
                thread_count++;
                ct << thread_count;

                //sub_bam = path_inter + "/thread_" + ct.str() + ".bam"; // write for each fasta file 
                ID_file = path_inter + "/IDs_" + ct.str() + ".txt";

                IDfile.open(ID_file.c_str());
                IDfile << *it_ID << endl; // IDs

                contig_counts = 2; // because first entry has already been added
                //ct.str(string()); //clear();
                
            }
        } else { // add remaining in last thread
            IDfile << *it_ID << endl; // IDs // ID start_pos end_pos

            contig_counts++;
        }
    }
    IDfile.close();
    // for last  ID file and last thread
    command_ROAST_mergeContigs_SCs = exe_path + "ROAST_mergeContigs_SCs " + ct.str() + " " + fastafile + " " + path_inter + " " + exe_path + " " + P1.str() + " " + P2.str() + " " + P3.str() + " " + P4.str() + " &";
    std::system(command_ROAST_mergeContigs_SCs.c_str());
    ct.str(string());


  /*  for (int i = 1; i <= allowed_threads; i++) {
        ct << i;
        string thread = ct.str();
       // ID_file = path_inter + "/IDs" + ct.str() + ".txt";

        //ROAST_mergeContigs_SCs_R2 tt;
        //tt.FindFragmentsBySCs(thread, fastafile,  path_inter, exe_path, SC_cons_len_TH, min_overlap_TH, blast_score_TH, max_allowed_gaps);
        command_ROAST_mergeContigs_SCs_R2 = exe_path + "ROAST_mergeContigs_SCs_R2 " + ct.str() + " " + fastafile + " " + path_inter + " " + exe_path + " " + P1.str() + " " + P2.str() + " " + P3.str() + " " + P4.str() + " &";
        std::system(command_ROAST_mergeContigs_SCs_R2.c_str());
        ct.str(string());
        //tms.FindFragmentsBySCs( ct.str(), extendfasta, path_inter);
    }*/

    P1.str(string());
    P2.str(string());
    P3.str(string());
    P4.str(string());
    running_threads_count = allowed_threads;

    while (running_threads_count > 0) {
        for (int i = 1; i <= allowed_threads; i++) {

            ct << i;
            string thread_done = path_inter + "/done_" + ct.str() + ".txt";
            ct.str(string());

            if (boost::filesystem::exists(thread_done)) // does filePath actually exist?
            {
                //continue;             
                running_threads_count--;

            } else { // unless exist status for all files obtained
                // i = 1; // start from 1st file again to check status
                sleep(10); // 
                running_threads_count = allowed_threads;
                break;
            }

        }
    }

   for (int i = 1; i <= allowed_threads; i++) {
        ct << i;
        string thread_done = path_inter + "/done_" + ct.str() + ".txt";
        utils.remove_file(thread_done.c_str());
        ct.str(string());

    }

    bySoftclip bsc;
    int merged_sc = bsc.broken_by_blast(new_assembly, iteration_log, fastafile);
    
    return merged_sc;
   //utils.remove_file(extendfasta);

}


int bySoftclip::merge_fragments(map<string, string> &final_merged_frags, map<string, string> &merged_fragments, string contig_1, int contig1_length, int contig2_length, string contig_2, bool extended_side, string updated_seq, string contig2_seq, int contig1_st, int contig1_end, int contig2_st, int contig2_end, bool strand_contig2, string header_fasta, string idd2, ofstream &log, set<string> &merged_contigs2, map<string, bool> &ID_SCside, bool &fragments_merged) {
    //contig1 is already merged
    
    // utils utils;
    typedef map<string, string> merged_fragment; // entries of this vectors are dir(read contig), updated id, updated seq, also separate for mate_contig
    typedef map<string, vector <string> > id_mapp;
    typedef map<string, string> final_merged_frag;
    bool extend_side;

    multimap<int, fragments_data>::reverse_iterator it_SC;
    string new_seq, merge, new_header_fasta;
    final_merged_frag::iterator it_update_seq;
    it_update_seq = final_merged_frags.find(header_fasta);

    int merged_count = 0;
    stringstream ss;

    if (!extended_side){ //RSC already extended
      // map contig1 SCs hit start and end positions on new updated contig for RSC / for LSC there is no change in st and end positions
        contig1_end = contig1_end + (updated_seq.size() - contig1_length);  // hypothetical 530 + (750-550) = 530 + 200 = 730
         contig1_st = contig1_end + (updated_seq.size() - contig1_length);  //              510 + (750-550) = 510 + 200 = 710 
        
    }
    if (!extended_side) { // RSC already extended so extended LSC now

        if (!strand_contig2) { // +
            // cout << "LSC +" << endl;

            merge = overlapmerge_left(contig2_seq, updated_seq, contig2_st, contig2_end, contig1_st);

            //merge_left = merge;

            if (merge != "") {

                ss << idd2 << "_and_" << header_fasta; // add [m] at the end or in new_header
                merged_count++;
                new_header_fasta = ss.str();

                if (new_header_fasta.length() > max_header_size) //27dec 21  to avoid error in cufflink due to long ID size
                    new_header_fasta = idd2;

                //28oct21 overlap boundry ectraction to avoid split of same position in multigene chimera
                int start_overlap = contig2_st; // start of overlap
                int end_overlap = contig2_end; // end of overlap
                merged_contigs.insert(make_pair(new_header_fasta, start_overlap));
                merged_contigs.insert(make_pair(new_header_fasta, end_overlap));


                merged_contigs2.insert(contig_2);
                ss.str(string());
                                log << contig_1 << "\t" << contig_2 << "\t" << "5' side" << "\t " << "+strand" << "\t" << start_overlap << "-" << end_overlap << endl;

            }
        } else { // if (strand_contig2) { // - / take reverse complement
            //   cout << " LSC -" << endl;

            string contig2_seq_RC = utils.Rcomplement(contig2_seq);
            int contig2_st_rc = contig2_length - contig2_end;
            int contig2_end_rc = contig2_length - contig2_st;
            
            merge = overlapmerge_left(contig2_seq_RC, updated_seq, contig2_st_rc, contig2_end_rc, contig1_st); // contig_st = 0
            //merge_left = merge;

            if (merge != "") {

                ss << idd2 << "_and_" << header_fasta; // add [m] at the end or in new_header
                merged_count++;
                new_header_fasta = ss.str();

                if (new_header_fasta.length() > max_header_size) //27dec 21  to avoid error in cufflink due to long ID size
                    new_header_fasta = idd2;

                ss.str(string());

                //28oct21 overlap boundry ectraction to avoid split of same position in multigene chimera
                int start_overlap = contig2_st_rc; // start of overlap
                int end_overlap = contig2_end_rc; // end of overlap
                merged_contigs.insert(make_pair(new_header_fasta, start_overlap));
                merged_contigs.insert(make_pair(new_header_fasta, end_overlap));
                
                // merged_contigs2.push_back(contig_1);
                merged_contigs2.insert(contig_2);

                log << contig_1 << "\t" << contig_2 << "\t" << "5' side" << "\t " << "-strand" << "\t" << start_overlap << "-" << end_overlap << endl;
                
            }

        }
        extend_side = true;
    } else { // LSC already extended so extended RSC now

        if (!strand_contig2) { // + 
            merge = overlapmerge_right(updated_seq, contig2_seq, contig1_end, contig2_st, contig2_end);
            
            if (merge != "") {

                ss << header_fasta << "_and_" << idd2; // << "_" << tool_name_flag << "_M_" << tool_name_flag; // "_merged"; // repeating things here, fixed*** // 
                merged_count++;
                new_header_fasta = ss.str();

                if (new_header_fasta.length() > max_header_size) //27dec 21  to avoid error in cufflink due to long ID size
                    new_header_fasta = header_fasta;

                //28oct21 overlap boundary extraction to avoid split of same position in multigene chimera
                int start_overlap = contig1_st; // start of overlap
                int end_overlap = contig1_end; // end of overlap
                merged_contigs.insert(make_pair(new_header_fasta, start_overlap));
                merged_contigs.insert(make_pair(new_header_fasta, end_overlap));
                
                //mergedAssembly << ">" << new_header << endl << merge << endl;
                log << contig_1 << "\t" << contig_2 << "\t" << "3' side" << "\t " << "+strand" << "\t" << start_overlap << "-" << end_overlap << endl;

                ss.str(string());
            }

        } else { // if (strand_contig2) { //-
            //  cout << "RSC - " << endl;
            string contig2_seq_RC = utils.Rcomplement(contig2_seq);
            int contig2_st_rc = contig2_length - contig2_end;
            int contig2_end_rc = contig2_length - contig2_st;
            
            merge = overlapmerge_right(updated_seq, contig2_seq_RC, contig1_end, contig2_st_rc, contig2_end_rc - 1);

            if (merge != "") {

                ss << header_fasta << "_and_" << idd2; // << "_" << tool_name_flag << "_M_" << tool_name_flag;
                merged_count++;
                new_header_fasta = ss.str();

                if (new_header_fasta.length() > max_header_size) //27dec 21  to avoid error in cufflink due to long ID size
                    new_header_fasta = header_fasta;

                //28oct21 overlap boundary extraction to avoid split of same position in multigene chimera
                int start_overlap = contig1_st; // start of overlap
                int end_overlap = contig1_end; // end of overlap
                merged_contigs.insert(make_pair(new_header_fasta, start_overlap));
                merged_contigs.insert(make_pair(new_header_fasta, end_overlap));
                
                ss.str(string());
                merged_contigs2.insert(contig_1);
                merged_contigs2.insert(contig_2);
                
                                // mergedAssembly << ">" << new_header << endl << merge << endl;
                log << contig_1 << "\t" << contig_2 << "\t" << "3' side" << "\t " << "-strand" << "\t" << start_overlap << "-" << end_overlap << endl;
            }
        }
        extend_side = false;
    }

    /*******updating old values********/
    // cout << "before updating values" << endl;
    if (!merge.empty()) {
        merged_fragments.insert(make_pair(contig_2, new_header_fasta)); // add new entry of contig_MI into merged_fragments

        //update id map with new id and add new entry in vector and erase the old one 
        ID_SCside.insert(make_pair(new_header_fasta, extend_side));
        //cout << "updated and erased id_map" << endl;

        // update final_merged_frags both new ids and updated seq and erase the old one

        final_merged_frags.insert(make_pair(new_header_fasta, merge));
        final_merged_frags.erase(it_update_seq);
        //cout << "updated and erased final_merged_frags" << endl;

        merged_contigs2.insert(contig_1);
        merged_contigs2.insert(contig_2); // to check for writing in final fasta file along with other unmerged contigs
        fragments_merged = true;
        /*******updating old values********/
    }

    return merged_count;

}

string bySoftclip::overlapmerge_left(string str1, string str2, int str1_st, int str1_end, int str2_st) {

    string final_seq = "";

    if ((str1.size() - str1_end) <= discard_contig_cor_len && str2_st <= discard_contig_cor_len) {

        final_seq = str1.substr(0, str1_st) + str2.substr(str2_st, str2.size() - 1); // left SC so overlapped region should be considered from str2(which is parent contig of SC)
    }
    else
    //    cout << str1.size() << " " << str1.size() - str1_end << "  " << str2.size() << " " << str2_st << endl;
    
    return final_seq;
}

string bySoftclip::overlapmerge_right(string str1, string str2, int str1_end, int str2_st, int str2_end) { //(contig1_seq, contig2_seq, contig1_end, contig2_st, contig2_end);


    string final_seq = "";

    if ((str1.size() - str1_end) <= discard_contig_cor_len && str2_st <= discard_contig_cor_len) {

        final_seq = str1.substr(0, str1_end) + str2.substr(str2_end, str2.size() - 1);
    }
    else
        //cout << str1.size() << " " << str1.size() - str1_end << "  " << str2.size() << " " << str2_st << endl;
    
    return final_seq;
}

/*31august 2021: same what BLAST is already doing, plus this function has error for while(j)*/
/*string bySoftclip::overlapmerge_left(string str1, string str2, int str1_st, int str2_st) {

    int discard_seq_TH = 10;
    int n, m, pos_MI = 1, pos_RI;
    int k = 1, consecutive_mm = 0;
    int counter = 0;
    int max_mismatch = 1, mismatch = 0;
    string overlap = "", maxOverlap = "", final, MI, RI, truncated_RI, truncated_MI;
    int overlapSize = 0;
    int i = str1_st, j = str2_st;
    char base_str1, base_str2;

    while (i < str1.size()) {

        while (j < str2.size()) {
            base_str1 = str1[i];
            base_str2 = str2[j];
            if (base_str1 == base_str2) {
                overlap = overlap + base_str1;
                i++;
                j++;
                consecutive_mm = 0;
                goto out_of_j;
            } else if ((j < str2.length() - 4) && (i < str1.length() - 4) && (str1[i] == str2[j + 1]) && (str1[i + 1] == str2[j + 2]) && (str1[i + 2] == str2[j + 3])) { // indel check, 3 bases check for conformation
                overlap = overlap + base_str2; //insertion in string2 
                j++;
                goto out_of_j;
            } else if ((j < str2.length() - 4) && (i < str1.length() - 4) && (str1[i + 1] == str2[j]) && (str1[i + 2] == str2[j + 1]) && (str1[i + 3] == str2[j + 2])) { // indel check
                overlap = overlap + base_str1; //insertion in string1 
                i++;
                goto out_of_j;
                pos_MI = j;
            } else {

                char code = base_code(base_str1, base_str2);
                mismatch++;
                if (mismatch <= max_mismatch && consecutive_mm <= 2) {
                    overlap = overlap + code;
                    consecutive_mm++;
                    i++;
                    j++;
                    pos_MI = j;

                    goto out_of_j;
                } else {
                    maxOverlap = overlap;
                    pos_MI = j; // size of overlap
                    pos_RI = i;
                    RI = str1.substr(0, str1_st);
                    MI = str2.substr(maxOverlap.size(), str2.size() - 1);
                    truncated_RI = str1.substr(pos_RI + 1, str1.size() - 1); //end of contig1 left from overlap
                    truncated_MI = str2.substr(pos_MI + 1, str2.size() - 1);
                    // std::transform(maxOverlap.begin(), maxOverlap.end(), maxOverlap.begin(), ::tolower); //end of contig2 left from overlap
                    final = RI + maxOverlap + MI;
                    overlap.clear();
                    goto out_of_i;
                }
            }
        }
out_of_j:
        ;
    }
out_of_i:
    ;
    if (!(overlap.empty())) { //when overlap found till end
        maxOverlap = overlap;
        RI = str1.substr(0, str1_st);
        MI = str2.substr(maxOverlap.size(), str2.size() - 1);
        truncated_RI = ""; //end of contig1 left from overlap
        truncated_MI = "";
        //std::transform(maxOverlap.begin(), maxOverlap.end(), maxOverlap.begin(), ::tolower); //end of contig2 left from overlap
        if (maxOverlap.size() < min_overlap_TH) {
            final = "";
            goto out;
        } else {
            final = RI + maxOverlap + MI;
            overlap.clear();
        }
    }
out:
    ;
    if (final != "") {

        if (maxOverlap.size() < min_overlap_TH) final = "";
        else {
            if ((truncated_RI.size() < discard_seq_TH) && (truncated_MI.size() <= truncated_RI.size())) {
                final = final +truncated_MI;
                // cout << pos_RI << endl << pos_MI << endl;
                //cout << "truncated RI:" << truncated_RI << endl;
                //cout << "truncated MI:" << truncated_MI << endl;
                return final;
            }//        else if ((truncated_MI.size() < discard_seq_TH) && (truncated_MI.size() < truncated_RI.size())) {

            else {
                //final = "";
                //cout << "no overlap found because truncated sequences exceeded specified range" << endl;
                return final;
            }
        }
    } else {
        // cout << "no overlap found" << endl;
        return final;
    }
}

string bySoftclip::overlapmerge_right(string str1, string str2, int str1_end, int str2_end) {

    int discard_seq_TH = 10;
    int pos_MI = 1, pos_RI;

    int counter = 0;
    int max_mismatch = 1, mismatch = 0, consecutive_mm = 0;
    string overlap = "", maxOverlap = "", final, MI, RI, truncated_RI, truncated_MI;
    int overlapSize = 0;
    int i = str1_end, j = str2_end - 1; //  start from 0
    char base_str1, base_str2;
    while (i >= 0) {
        while (j >= 0) {
            if (j == 0) goto out_of_i;
            base_str1 = str1[i];
            base_str2 = str2[j];
            //    for (k = 1; k <= str1.size() && k <= str2.size(); k++)
            //   {
            if (base_str1 == base_str2) {
                overlap = base_str1 + overlap;
                i--;
                j--;
                consecutive_mm = 0;
                goto out_of_j;
            } else if ((j > 3) && (i > 3) && (str1[i] == str2[j - 1]) && (str1[i - 1] == str2[j - 2]) && (str1[i - 2] == str2[j - 3])) { // indel check, 3 bases check for conformation
                overlap = base_str2 + overlap; //insertion in string2 
                j--;
                pos_MI = j;
                goto out_of_j;
            } else if ((j < str2.length() - 4) && (i < str1.length() - 4) && (str1[i - 1] == str2[j]) && (str1[i - 2] == str2[j - 1]) && (str1[i - 3] == str2[j - 2])) { // indel check
                overlap = base_str1 + overlap; //insertion in string1 
                i--;
                pos_MI = j;
                goto out_of_j;
            } else {

                char code = base_code(base_str1, base_str2);
                mismatch++;
                if (mismatch <= max_mismatch && consecutive_mm <= 2) {
                    overlap = code + overlap;
                    consecutive_mm++;
                    i--;
                    j--;
                    pos_MI = j;
                    goto out_of_j;
                } else {
                    maxOverlap = overlap;
                    pos_MI = j; // size of overlap
                    pos_RI = str1.size() - overlap.size();
                    RI = str1.substr(0, pos_RI);
                    MI = str2.substr(str2_end, str2.size() - 1);
                    truncated_RI = ""; //str1.substr(0, pos_RI +); //end of contig1 left from overlap
                    truncated_MI = str2.substr(0, pos_MI);
                    //std::transform(maxOverlap.begin(), maxOverlap.end(), maxOverlap.begin(), ::tolower); //end of contig2 left from overlap
                    final = RI + maxOverlap + MI;
                    overlap.clear();
                    goto out_of_i;
                }
            }
        }

out_of_j:
        ;
    }
out_of_i:
    ;
    if (!(overlap.empty())) { //when overlap found till end
        maxOverlap = overlap;
        pos_RI = str1.size() - overlap.size();
        RI = str1.substr(0, pos_RI - 1); //contig 1 from start till start of overlap
        MI = str2.substr(str2_end, str2.size() - 1); //contig2 from end of overlap till end of contig
        truncated_RI = "";
        truncated_MI = "";
        std::transform(maxOverlap.begin(), maxOverlap.end(), maxOverlap.begin(), ::tolower); //end of contig2 left from overlap

        if (maxOverlap.size() < min_overlap_TH) {
            final = "";
            goto out;
        } else {
            final = RI + maxOverlap + MI;
            overlap.clear();
        }
        // goto out_of_i;
    }
out:
    ;
    if (final != "") {

        if (maxOverlap.size() < min_overlap_TH) final = "";
        else {
            if ((truncated_RI.size() < discard_seq_TH) && (truncated_MI.size() <= truncated_RI.size())) {//if ((truncated_RI.size() < discard_seq_TH)) {
                // final = final + truncated_MI;
                //cout << "truncated RI:" << truncated_RI << endl;
                // cout << "truncated MI:" << truncated_MI << endl;
                return final;
            }//        else if ((truncated_MI.size() < discard_seq_TH) && (truncated_MI.size() < truncated_RI.size() )) {

            else {
                final = "";
                // cout << "no overlap found beacuse truancted sequences exceeded specified range" << endl;
                return final;
            }
        }
    } else {
        //cout << "no overlap found" << endl;
        return final;
    }
} */

char bySoftclip::base_code(char base1, char base2) {
    char code;
    if ((base1 == 'G' && base2 == 'A') || (base1 == 'A' && base2 == 'G')) code = 'R'; // puRine
    else if ((base1 == 'T' && base2 == 'C') || (base1 == 'C' && base2 == 'T')) code = 'Y'; //pYrimidine
    else if ((base1 == 'A' && base2 == 'C') || (base1 == 'C' && base2 == 'A')) code = 'M'; //aMino
    else if ((base1 == 'G' && base2 == 'T') || (base1 == 'T' && base2 == 'G')) code = 'K'; //Keto
    else if ((base1 == 'G' && base2 == 'C') || (base1 == 'C' && base2 == 'G')) code = 'S'; //Strong interaction (3 H bonds)
    else if ((base1 == 'A' && base2 == 'T') || (base1 == 'T' && base2 == 'A')) code = 'W'; //Weak interaction (2 H bonds)
    else code = 'N';
    return code;


}
