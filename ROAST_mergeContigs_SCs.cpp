
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ROAST_mergeContig_SCs.cpp
 * Author: madiha
 * 
 * Created on March 31, 2021, 2:56 PM
 */

#include "ROAST_mergeContigs_SCs.h"
#include "utils.h"
#include "global.h"
#include "alignment.h"
#include "bySoftclip.h"
#include <sstream>
#include <fstream>
#include <iostream>                                                                                                                                                                                          
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
#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>

using namespace std;

ROAST_mergeContigs_SCs::ROAST_mergeContigs_SCs() {
}

ROAST_mergeContigs_SCs::ROAST_mergeContigs_SCs(const ROAST_mergeContigs_SCs& orig) {
}

ROAST_mergeContigs_SCs::~ROAST_mergeContigs_SCs() {
}

//string delimiter = "\t";

void ROAST_mergeContigs_SCs::FindFragmentsBySCs(string thread, string extended_assembly, string path_inter, string exe_path, int SC_cons_len_TH, int min_overlap_TH, int blast_score_TH, int max_allowed_gaps) {
    string delimiter = "\t";

    utils utils;

    string ID_file = path_inter + "/IDs_" + thread + ".txt";

    string filterd_frags = path_inter + "/filtered_frags_" + thread + ".txt";

    string SC_query_left = path_inter + "/SC_query_left_" + thread;
    string SC_query_right = path_inter + "/SC_query_right_" + thread;

    string SC_subject_left = path_inter + "/SC_sub_left_" + thread;
    string SC_subject_right = path_inter + "/SC_sub_right_" + thread;

    string blast_output_left = path_inter + "/blast_output_left_" + thread;
    string blast_output_right = path_inter + "/blast_output_right_" + thread;

    string blast_output_left2 = path_inter + "/blast_output_left2_" + thread;
    string blast_output_right2 = path_inter + "/blast_output_right2_" + thread;



    string strand, strand2, SC_seq, contig1_seq, contig2_seq, contig2_seq_RC, ID, merge, id1, id2, query_seq;
    int contig2_st, contig1_st, contig2_end, contig1_end, contig1_length, contig2_length, contig2_st_rc, contig2_end_rc;
    string idd1, idd2, contig_1, contig_2;

    bool SC_side, consider_fragment, check_BLAST = true;
    vector<string> temp, same_contig;

    string entry2, contig_name, contig_left_merged, entry;
    int count = 0, merged_count = 0, query_size_check = min_overlap_TH + 5;

    set<string> IDs;

    ofstream filt_frags;
    filt_frags.open(filterd_frags.c_str());

    ifstream in_IDs;
    in_IDs.open(ID_file.c_str());

    while (!getline(in_IDs, ID).eof()) { //save all IDs of fastafile
        IDs.insert(ID);
    }
    in_IDs.close();


    for (set<string>::iterator it_ID = IDs.begin(); it_ID != IDs.end(); it_ID++) {

        string contig_name = *it_ID;

        if (contig_name.find("trLE") != string::npos || contig_name.find("trBE") != string::npos) { // consider if Left side extended or both


            ofstream fasta_query_left, fasta_subject_left;
            //cout << "run left blast" << endl;
            query_seq = utils.extract_fasta(*it_ID, extended_assembly);

            if (query_seq.size() > min_overlap_TH + 5) {

                fasta_query_left.open(SC_query_left.c_str());
                fasta_query_left << ">" << *it_ID << endl << query_seq.substr(0, min_overlap_TH) << endl; // for those extended contigs which extended >= min_overlap_TH in multiple iteration got check here at the end of inner iterations
                fasta_query_left.close();

                //blastn -evalue 5 -word_size 4 -outfmt '7 qseqid qlen sseqid slen qcovs frames pident nident length mismatch qstart qend sstart send' -subject -query  -out
                string blast_command_left = "blastn -evalue 1 -word_size 4 -outfmt '7 qseqid qlen sseqid slen qcovhsp nident length mismatch gaps qstart qend sstart send sstrand' -subject " + extended_assembly + " -query " + SC_query_left + " -out " + blast_output_left + "  > /dev/null 2>&1";
                std::system(blast_command_left.c_str());

                ifstream blast_left(blast_output_left.c_str());

                while (!(getline(blast_left, entry).eof())) {
                    //cout << "left BLAT" << endl;
                    if (entry[0] == '#') {
                        continue;
                    } else {
                        temp.clear();
                        utils.str_split(entry, temp, delimiter);

                        string id1 = temp[0];
                        string id2 = temp[2];
                        float query_size = atoi(temp[1].c_str());
                        //int target_length = abs(target_end - target_start);
                        float score1 = atoi(temp[5].c_str()); //nident
                        int gaps = atoi(temp[8].c_str());

                        int final_score = float(score1 / query_size) * 100; // score1 is number of matches

                        if (id1 == id2) // || g_id_read == g_id_mate)) {
                            consider_fragment = false;

                        else if (id1.find(id2) != string::npos && id2.find(id1) != string::npos) //read contig is not part of mate contig and vice versa *logically it shouldn't be the case so check that*
                            consider_fragment = false;

                        else
                            consider_fragment = true;
                        //first condition because query soft clipped is already stitched in its own sequence

                        if (consider_fragment && gaps <= max_allowed_gaps && final_score >= blast_score_TH) {// && score1 >= min_overlap_TH) {

                            int target_size = atoi(temp[3].c_str());
                            int target_start = atoi(temp[11].c_str());
                            int target_end = atoi(temp[12].c_str());

                            contig1_st = atoi(temp[9].c_str());
                            contig1_end = atoi(temp[10].c_str());


                            if (temp[13] == "plus")
                                strand = "+";
                            else // target_start > target_end
                                strand = "-";

                            stringstream entry_new;
                            /*second BLAST between left-over corner of contigB (if>= 10bases) and remaining sequence of contigA
                             1. for < 10 bases consider them fragment and first BLAST hit as final
                             2. for >= 10 bases if hit starts from boundary of contigA consider them as fragment
                             3. else don't consider them fragment and ignore */


                            string id2_seq = utils.extract_fasta(id2, extended_assembly);
                            string rem_bases_id2, rem_bases_id1, query2_seq;

                            if (strand == "+") {
                                rem_bases_id2 = id2_seq.substr(target_end, id2_seq.size() - 1);
                                query2_seq = id2_seq.substr(target_start, id2_seq.size() - 1); // take whole fragment of Target1 again for BLAST2
                            } else {
                                rem_bases_id2 = id2_seq.substr(0, target_end);
                                query2_seq = id2_seq.substr(0, target_start); // take whole fragment of Target1 again for BLAST2
                            }

                            if (rem_bases_id2.size() >= 10) { // BLAST2 only if 

                                fasta_subject_left.open(SC_subject_left.c_str());
                                fasta_subject_left << ">" << id1 << endl << query_seq << endl; // if fragment found from first BLAST hit do BLAST2 for full contigA and whole fragment of contigB
                                fasta_subject_left.close();

                                fasta_query_left.open(SC_query_left.c_str());
                                fasta_query_left << ">" << id2 << endl << query2_seq << endl; // if fragment found from first BLAST hit do BLAST2 for remaining portion of contigA and contigB
                                fasta_query_left.close();

                                string blast_command_left = "blastn -evalue 1 -word_size 4 -outfmt '7 qseqid qlen sseqid slen qcovhsp nident length mismatch gaps qstart qend sstart send sstrand' -subject " + SC_subject_left + " -query " + SC_query_left + " -out " + blast_output_left2 + "  > /dev/null 2>&1";
                                std::system(blast_command_left.c_str());

                                ifstream blast_left2(blast_output_left2.c_str());

                                while (!(getline(blast_left2, entry2).eof())) {
                                    /*NOw target is contig A and query is contig B
                                         update position of contigA(target-> previous query) as contig B (query -> previous target) is already till corner 
                                     */
                                    //cout << "left BLAT" << endl;
                                    if (entry2[0] == '#') {
                                        check_BLAST = true;
                                        continue;
                                    } else if (check_BLAST) {
                                        temp.clear();
                                        utils.str_split(entry2, temp, delimiter);
                                        
                                        if (temp[13] == "plus")
                                            strand2 = "+";
                                        else // target_start > target_end
                                            strand2 = "-";

                                        if (strand == strand2) { // 
                                                //string id1 = temp[0];
                                                //string id2 = temp[2];
                                                float query_size = atoi(temp[1].c_str());
                                                //int target_length = abs(target_end - target_start);
                                                float score1 = atoi(temp[5].c_str()); //nident
                                                int gaps = atoi(temp[8].c_str());

                                            int final_score = float(score1 / query_size) * 100; // score1 is number of matches

                                            /*double check on query size and score TH because of 190 FP in rice first iteration*/
                                            if (query_size >= query_size_check) {
                                                //blast_score_TH = 90; // keep it same
                                            } else if (query_size >= min_overlap_TH) {
                                                blast_score_TH = 100;
                                            }
                                            if (gaps <= max_allowed_gaps && final_score >= blast_score_TH) {
                                                // int target_size = atoi(temp[3].c_str());
                                                int contigA_target_st = atoi(temp[11].c_str());
                                                int contigA_target_en = atoi(temp[12].c_str());

                                                int contigB_query_st = atoi(temp[9].c_str());
                                                int contigB_query_end = atoi(temp[10].c_str());

                                                if (strand2 == "+") {

                                                    contig1_st = contigA_target_st;
                                                    contig1_end = contigA_target_en;
                                                    //target_start = would be same as found in BLAST1
                                                    target_end = target_start + contigB_query_end; // 450 + 45 = 495 (length 500)
                                                    // strand = "+"; use strand of BLAST1, don't update here 

                                                    rem_bases_id2 = id2_seq.substr(target_end, id2_seq.size() - 1);
                                                    rem_bases_id1 = query_seq.substr(0, contig1_st);

                                                } else { // to make Q1(target2) as positive and Q2(target1) as negative as in BLAST1

                                                    contig1_st = query_seq.size() - contigA_target_en;
                                                    contig1_end = query_seq.size() - contigA_target_st;
                                                    target_start = contigB_query_end; // 30 for length 30
                                                    target_end = contigB_query_st; // 5 for length of 30
                                                    //strand = "-";

                                                    rem_bases_id2 = id2_seq.substr(0, target_end);
                                                    rem_bases_id1 = query_seq.substr(0, contig1_st);
                                                }

                                                //first condition because query soft clipped is already stitched in its own sequence
                                                if (rem_bases_id1.size() < 10 && rem_bases_id2.size() < 10) {// && score1 >= min_overlap_TH) {

                                                    stringstream entry_new;
                                                    entry_new << "left" << "\t" << id1 << "\t" << id2 << "\t" << final_score << "\t" << strand << "\t" << contig1_st << "\t" << contig1_end << "\t" << target_size << "\t" << target_start << "\t" << target_end << "\t" << "0" << "\t" << query_size;
                                                    string blast_hits = entry_new.str();
                                                    same_contig.push_back(blast_hits);

                                                    //cout << entry2 << endl;
                                                    entry_new.str(string());
                                                } else
                                                    check_BLAST = false; // check only top hit of BLAST to reduce time
                                            }
                                        } // strand check.. should be same in first and second BLAST 14jan21
                                    }
                                    consider_fragment = false;
                                    break; // only first best hit
                                }
                                blast_left2.close();


                            } else if (final_score == 100) { //100// remaining bases are < discard_contig_cor_len so seed should have 100% match

                                entry_new << "left" << "\t" << id1 << "\t" << id2 << "\t" << final_score << "\t" << strand << "\t" << contig1_st << "\t" << contig1_end << "\t" << target_size << "\t" << target_start << "\t" << target_end << "\t" << "0" << "\t" << query_size;
                                string blast_hits = entry_new.str();
                                same_contig.push_back(blast_hits);
                                consider_fragment = false;
                                // break; //only first best hit
                            }
                            goto next_contig_hit_left; // to go to next subject hit 
                        }
                    }
next_contig_hit_left:
                    ;
                    consider_fragment = false;

                }
                blast_left.close();


                //cout << "before filter" << endl;
                /*add script from merge_contigs_sc.cpp 
                It filters hits for each SC blast and select best one to select for merging*/
                vector<string> contig2_list, same_contig_filtered, temp3, temp4;
                same_contig_filtered.reserve(10);
                contig2_list.reserve(10);


                for (size_t itr_same = 0; itr_same < same_contig.size(); itr_same++) { //separate contig2 list to filter repetition

                    string str = same_contig[itr_same];
                    vector<string> entry_split;

                    utils.str_split(str, entry_split, delimiter);
                    contig2_list.push_back(entry_split[2]); // collect mate contig
                }

                for (size_t itr_list = 0; itr_list < contig2_list.size(); itr_list++) {

                    string name = contig2_list[itr_list];
                    int mycount = std::count(contig2_list.begin(), contig2_list.end(), name); //filter for same mate contigs for left and right

                    if (mycount > 1) // erase that contig and dont add in new vector;
                        continue;
                    else //if not found add whole entry in new vector to process further               
                        same_contig_filtered.push_back(same_contig[itr_list]);

                } // filtering done

                float max_score_left = 0.0, max_score_right = 0.0;
                string max_temp_left, max_temp_right;
                float current_score = 0.0;
                stringstream ss;
                map<string, string>::iterator it_fasta;
                //      cout << "after filter" << endl;

                for (int itr = 0; itr < same_contig_filtered.size(); itr++) {

                    string temp2 = same_contig_filtered[itr];
                    temp3.clear();
                    utils.str_split(temp2, temp3, delimiter);

                    if (temp3[0] == "left")
                        SC_side = true; //LSC
                    else
                        SC_side = false; //RSC

                    id1 = temp3[1];
                    id2 = temp3[2];
                    current_score = atof(temp3[3].c_str());
                    strand = temp3[4];
                    contig2_st = atoi(temp3[8].c_str());
                    contig2_end = atoi(temp3[9].c_str()) - 1;
                    contig2_length = atoi(temp3[7].c_str());
                    //SC_seq = temp3[11];

                    //std::size_t contig1_st = contig1_seq.find(SC_seq);
                    contig1_seq = utils.extract_fasta(id1, extended_assembly);

                    contig1_end = contig1_seq.size() - 1;
                    contig1_length = contig1_seq.size();

                    int TH;
                    if (contig2_seq.size() <= 300)
                        TH = 100;
                    else TH =
                            200;
                    // Check for max score for left and right ignoring positions at wrong strand
                    //          cout << "before check max" << endl;
                    // if (SC_side) { //LSC
                    if (strand == "+") {

                        if ((contig2_end >= contig2_length - TH) && ((contig2_length - contig2_st) < contig1_length)) {
                            if (max_score_left < current_score) {
                                max_temp_left = temp2;
                                max_score_left = current_score;
                            }
                        }
                        //cout << id1 << "  " << id2 << endl;
                    } else {// if (strand == "-") {
                        contig2_seq_RC = utils.Rcomplement(contig2_seq);

                        contig2_st_rc = contig2_length - contig2_end;
                        contig2_end_rc = contig2_length - contig2_st;

                        if ((contig2_end_rc >= contig2_length - TH) && ((contig2_length - contig2_st_rc) < contig1_length)) {
                            if (max_score_left < current_score) {
                                max_temp_left = temp2;
                                max_score_left = current_score;
                            }
                        }
                        // cout << id1 << "  " << id2 << endl;
                    }
                    //  }

                }
                //   cout << "after check max" << endl;
                if (!(max_temp_left.empty())) {

                    filt_frags << max_temp_left << endl;
                    //broken_bySCs.insert(make_pair(fragments_data.BLAT_score, fragments_data));

                }
                //*********************************
                same_contig.clear();
                same_contig_filtered.clear();
            } // query size > min_overlap_TH
        }
        if (contig_name.find("trRE") != string::npos || contig_name.find("trBE") != string::npos) { // consider if right side is extended or both sides

            ofstream fasta_query_right, fasta_subject_right;
            //cout << "run left blast" << endl;
            query_seq = utils.extract_fasta(contig_name, extended_assembly);

            if (query_seq.size() > min_overlap_TH) { // 5/1/22 because of query size less than overlap TH

                fasta_query_right.open(SC_query_right.c_str());
                fasta_query_right << ">" << contig_name << endl << query_seq.substr(query_seq.size() - min_overlap_TH, query_seq.size()) << endl; // for those extended contigs which extended >= min_overlap_TH in multiple iteration got check here at the end of inner iterations
                fasta_query_right.close();
                // string blast_command_left = exe_path + "external_tools/blast -stepSize=5 -repMatch=2253 -minScore=40 -minIdentity=50 " + extended_assembly + " " + SC_query_left + " " + blast_output_left + "  > /dev/null 2>&1";
                // std::system(blast_command_left.c_str());

                //blastn -evalue 5 -word_size 4 -outfmt '7 qseqid qlen sseqid slen qcovs frames pident nident length mismatch qstart qend sstart send' -subject -query  -out
                string blast_command_right = "blastn -evalue 1 -word_size 4 -outfmt '7 qseqid qlen sseqid slen qcovhsp nident length mismatch gaps qstart qend sstart send sstrand' -subject " + extended_assembly + " -query " + SC_query_right + " -out " + blast_output_right + "  > /dev/null 2>&1";
                std::system(blast_command_right.c_str());

                ifstream blast_right(blast_output_right.c_str());

                while (!(getline(blast_right, entry).eof())) {
                    //cout << "left right" << endl;
                    if (entry[0] == '#') {
                        continue;
                    } else {
                        temp.clear();
                        utils.str_split(entry, temp, delimiter);

                        string id1 = temp[0];
                        string id2 = temp[2];
                        float query_size = atoi(temp[1].c_str());

                        //int target_length = abs(target_end - target_start);
                        float score1 = atoi(temp[5].c_str()); //nident
                        int gaps = atoi(temp[8].c_str());

                        int final_score = float(score1 / query_size) * 100; // score1 is number of matches

                        if (id1 == id2) // || g_id_read == g_id_mate)) {
                            consider_fragment = false;

                        else if (id1.find(id2) != string::npos && id2.find(id1) != string::npos) //read contig is not part of mate contig and vice versa *logically it shouldn't be the case so check that*
                            consider_fragment = false;

                        else
                            consider_fragment = true;

                        if (consider_fragment && gaps <= max_allowed_gaps && final_score >= blast_score_TH) {// && score1 >= min_overlap_TH) {

                            int target_size = atoi(temp[3].c_str());
                            int target_start = atoi(temp[11].c_str());
                            int target_end = atoi(temp[12].c_str());

                            contig1_st = atoi(temp[9].c_str());
                            contig1_end = atoi(temp[10].c_str());


                            if (temp[13] == "plus")
                                strand = "+";
                            else // target_start > target_end
                                strand = "-";


                            //first condition because query soft clipped is already stitched in its own sequence
                            stringstream entry_new;
                            /*second BLAST between left-over corner of contigB (if>= 10bases) and remaining sequence of contigA
                             1. for < 10 bases consider them fragment and first BLAST hit as final
                             2. for >= 10 bases if hit starts from boundary of contigA consider them as fragment
                             3. else don't consider them fragment and ignore */

                            string id2_seq = utils.extract_fasta(id2, extended_assembly);
                            string rem_bases_id2, rem_bases_id1, query_seq2;

                            if (strand == "+") {
                                rem_bases_id2 = id2_seq.substr(0, target_start);
                                query_seq2 = id2_seq.substr(0, target_end);

                            } else {
                                rem_bases_id2 = id2_seq.substr(target_start, id2_seq.size() - 1);
                                query_seq2 = id2_seq.substr(target_end, id2_seq.size() - 1);
                            }

                            if (rem_bases_id2.size() >= 10) { // BLAST2 only if 

                                fasta_subject_right.open(SC_subject_right.c_str());
                                fasta_subject_right << ">" << id1 << endl << query_seq << endl; // if fragment found from first BLAST hit do BLAST2 for whole of contigA and and full fragment of contigB
                                fasta_subject_right.close();

                                fasta_query_right.open(SC_query_right.c_str());
                                fasta_query_right << ">" << id2 << endl << query_seq2 << endl; // if fragment found from first BLAST hit do BLAST2 for remaining portion of contigA and contigB
                                fasta_query_right.close();

                                string blast_command_right = "blastn -evalue 1 -word_size 4 -outfmt '7 qseqid qlen sseqid slen qcovhsp nident length mismatch gaps qstart qend sstart send sstrand' -subject " + SC_subject_right + " -query " + SC_query_right + " -out " + blast_output_right2 + "  > /dev/null 2>&1";
                                std::system(blast_command_right.c_str());

                                ifstream blast_right2(blast_output_right2.c_str());

                                while (!(getline(blast_right2, entry2).eof())) {
                                    /*NOw target is contig A and query is contig B
                                     update position of contigA(target-> previous query) as contig B (query -> previous target) is already till corner 
                                     */
                                    //cout << "left BLAT" << endl;
                                    if (entry2[0] == '#') {
                                        check_BLAST = true;
                                        continue;
                                    } else if (check_BLAST) {
                                        temp.clear();
                                        utils.str_split(entry2, temp, delimiter);

                                        //string id1 = temp[0];
                                        //string id2 = temp[2];
                                        float query_size = atoi(temp[1].c_str());

                                        if (temp[13] == "plus")
                                            strand2 = "+";
                                        else // target_start > target_end
                                            strand2 = "-";

                                        if (strand == strand2) { //

                                            //int target_length = abs(target_end - target_start);
                                            float score1 = atoi(temp[5].c_str()); //nident
                                            int gaps = atoi(temp[8].c_str());

                                            /*double check on query size and score TH because of 190 FP in rice first iteration*/
                                            if (query_size >= query_size_check) {
                                                //blast_score_TH = 90; // keep it same
                                            } else if (query_size >= min_overlap_TH + 5) {
                                                blast_score_TH = 100;
                                            }

                                            int final_score = float(score1 / query_size) * 100; // score1 is number of matches

                                            if (gaps <= max_allowed_gaps && final_score >= blast_score_TH) {
                                                // int target_size = atoi(temp[3].c_str());
                                                int contigA_target_st = atoi(temp[11].c_str());
                                                int contigA_target_en = atoi(temp[12].c_str());

                                                int contigB_query_st = atoi(temp[9].c_str());
                                                int contigB_query_end = atoi(temp[10].c_str());

                                                if (strand2 == "+") {

                                                    contig1_st = contigA_target_st;
                                                    contig1_end = contigA_target_en; //contig1_end = atoi(temp[10].c_str()); // 20 + 40 (target2 end)
                                                    target_start = contigB_query_st;
                                                    target_end = contigB_query_end;
                                                    // strand = "+"; use strand of BLAST1, don't update here 

                                                    rem_bases_id2 = id2_seq.substr(0, target_start);
                                                    rem_bases_id1 = query_seq.substr(contig1_end, query_seq.size() - 1);

                                                } else { // target_start > target_end

                                                    contig1_st = query_seq.size() - contigA_target_en;
                                                    contig1_end = query_seq.size() - contigA_target_st;
                                                    target_start = target_end + contigB_query_end;
                                                    // target_end same as in BLAST1 
                                                    //strand = "-";

                                                    rem_bases_id2 = id2_seq.substr(target_start, id2_seq.size() - 1);
                                                    rem_bases_id1 = query_seq.substr(contig1_end, query_seq.size() - 1);
                                                }

                                                //first condition because query soft clipped is already stitched in its own sequence
                                                if (rem_bases_id1.size() < 10 && rem_bases_id2.size() < 10) {// && score1 >= min_overlap_TH) {

                                                    stringstream entry_new;
                                                    entry_new << "right" << "\t" << id1 << "\t" << id2 << "\t" << final_score << "\t" << strand << "\t" << contig1_st << "\t" << contig1_end << "\t" << target_size << "\t" << target_start << "\t" << target_end << "\t" << "1" << "\t" << query_size;
                                                    string blast_hits = entry_new.str();
                                                    same_contig.push_back(blast_hits);

                                                    //cout << entry2 << endl;
                                                    entry_new.str(string());
                                                    break; // only first best hit

                                                } else
                                                    check_BLAST = false; // check only top entry to reduce time 14jan21
                                            }
                                        } // strand check.. should be same in first and second BLAST 14jan21

                                    }
                                    consider_fragment = false;

                                }
                                blast_right2.close();

                                // break; //  consider only first hit
                            } else if (final_score == 100) { //100 // remaining bases are < discard_contig_cor_len // if overlap is of seed size get 100 % match

                                entry_new << "right" << "\t" << id1 << "\t" << id2 << "\t" << final_score << "\t" << strand << "\t" << contig1_st << "\t" << contig1_end << "\t" << target_size << "\t" << target_start << "\t" << target_end << "\t" << "1" << "\t" << query_size;
                                string blast_hits = entry_new.str();
                                same_contig.push_back(blast_hits);
                                consider_fragment = false;
                                // break; // only first best hit
                            }
                            goto next_contig_hit_right;
                        }

                    }
next_contig_hit_right:
                    ;
                    consider_fragment = false;
                }
                blast_right.close();


                /*add script from merge_contigs_sc.cpp 
                       It filters hits for each SC blast and select best one to select for merging*/
                vector<string> contig2_list, same_contig_filtered, temp3, temp4;
                same_contig_filtered.reserve(10);
                contig2_list.reserve(10);

                //cout << "before filter" << endl;
                for (size_t itr_same = 0; itr_same < same_contig.size(); itr_same++) { //separate contig2 list to filter repetition

                    string str = same_contig[itr_same];
                    vector<string> entry_split;

                    utils.str_split(str, entry_split, delimiter);
                    contig2_list.push_back(entry_split[2]); // collect mate contig
                }

                for (size_t itr_list = 0; itr_list < contig2_list.size(); itr_list++) {

                    string name = contig2_list[itr_list];
                    int mycount = std::count(contig2_list.begin(), contig2_list.end(), name); //filter for same mate contigs for left and right

                    if (mycount > 1) // erase that contig and dont add in new vector;
                        continue;
                    else //if not found add whole entry in new vector to process further               
                        same_contig_filtered.push_back(same_contig[itr_list]);

                } // filtering done

                //      cout << "after filter" << endl;
                float max_score_left = 0.0, max_score_right = 0.0;
                string max_temp_left, max_temp_right;
                float current_score = 0.0;
                stringstream ss;
                map<string, string>::iterator it_fasta;
                for (int itr = 0; itr < same_contig_filtered.size(); itr++) {

                    string temp2 = same_contig_filtered[itr];
                    temp3.clear();
                    utils.str_split(temp2, temp3, delimiter);

                    if (temp3[0] == "left")
                        SC_side = true; //LSC
                    else
                        SC_side = false; //RSC

                    id1 = temp3[1];
                    id2 = temp3[2];
                    current_score = atof(temp3[3].c_str());
                    strand = temp3[4];
                    contig2_st = atoi(temp3[8].c_str());
                    contig2_end = atoi(temp3[9].c_str()) - 1;
                    contig2_length = atoi(temp3[7].c_str());
                    //SC_seq = temp3[11];

                    //std::size_t contig1_st = contig1_seq.find(SC_seq);
                    contig1_seq = utils.extract_fasta(id1, extended_assembly);

                    contig1_end = contig1_seq.size() - 1;
                    contig1_length = contig1_seq.size();

                    int TH;
                    if (contig2_seq.size() <= 300)
                        TH = 100;
                    else TH =
                            200;
                    // Check for max score for left and right ignoring positions at wrong strand
                    //          cout << "before check max" << endl;

                    // if (!SC_side) { //RSC
                    if (strand == "+") {
                        if ((contig2_st <= TH) && (contig2_end < contig1_length)) {
                            if (max_score_right < current_score) {
                                max_temp_right = temp2;
                                max_score_right = current_score;
                            }
                        }
                        //  cout << id1 << "  " << id2 << endl;
                    } else { //if (strand == "-") {

                        contig2_seq_RC = utils.Rcomplement(contig2_seq);
                        contig2_st_rc = contig2_length - contig2_end;
                        contig2_end_rc = contig2_length - contig2_st;

                        if ((contig2_st_rc <= TH) && (contig2_end_rc < contig1_length)) {

                            if (max_score_right < current_score) {
                                max_temp_right = temp2;
                                max_score_right = current_score;
                            }
                        }
                    }
                    //  }
                    //}

                    if (!(max_temp_right.empty())) {
                        filt_frags << max_temp_right << endl;
                        //broken_bySCs.insert(make_pair(fragments_data.BLAT_score, fragments_data));
                    }//*********************************
                    same_contig.clear();
                    same_contig_filtered.clear();
                }
            } // query size > overlap_TH
        }
    }
    filt_frags.close();
    utils.remove_file(SC_query_left.c_str());
    utils.remove_file(SC_query_right.c_str());

    utils.remove_file(blast_output_left.c_str());
    utils.remove_file(blast_output_right.c_str());

    utils.remove_file(SC_subject_left.c_str());
    utils.remove_file(SC_subject_right.c_str());

    utils.remove_file(blast_output_left2.c_str());
    utils.remove_file(blast_output_right2.c_str());
    utils.remove_file(ID_file.c_str());

    string thread_done = path_inter + "/done_" + thread + ".txt";
    ofstream done;
    done.open(thread_done.c_str());
    done << "Job Finished";
    done.close();

    /* string new_assembly = path_inter + "/Improvrd_assembly_SCsALL.fasta";
     string iteration_log = path_inter + "/Improvrd_assembly_SCsALL.log";

     bySoftclip bsc;
     bsc.broken_by_blast(new_assembly, iteration_log, extended_assembly, merged_contigs);*/
}

int main(int argc, char** argv) {

    if (argc < 2) {
        cout << "No arguments found" << endl;
        exit(0);
    }
    vector <string> args(argv, argv + argc);

    string thread = args[1];

    string extendedfasta = args[2];

    string path_inter = args[3];

    int SC_cons_len_TH = atoi(args[4].c_str());

    string exe_path = args[5];

    int min_overlap_TH = atoi(args[6].c_str());

    int blast_score_TH = atoi(args[7].c_str());

    int max_allowed_gaps = atoi(args[8].c_str());

    ROAST_mergeContigs_SCs tt;
    tt.FindFragmentsBySCs(thread, extendedfasta, path_inter, exe_path, SC_cons_len_TH, min_overlap_TH, blast_score_TH, max_allowed_gaps);


    //cout << "ROAST_mergeContigs_SCs_R2 returned successfully " << endl;

    return 0;
}
