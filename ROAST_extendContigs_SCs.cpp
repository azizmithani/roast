/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ROAST_extendContigs_SCs.cpp
 * Author: madiha
 * 
 * Created on January 7, 2021, 12:24 PM
 */

#include "ROAST_extendContigs_SCs.h"
#include "utils.h"
#include <sstream>
#include <iostream>
#include <fstream>                                                                                                                                                                                         
#include <vector>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <cctype>    
#include <climits>
#include <set> 
#include <map>
#include <list>
#include <api/BamAlignment.h>
#include <api/BamReader.h>
using namespace BamTools;
using namespace std;

ROAST_extendContigs_SCs::ROAST_extendContigs_SCs() {
}

ROAST_extendContigs_SCs::ROAST_extendContigs_SCs(const ROAST_extendContigs_SCs& orig) {
}

ROAST_extendContigs_SCs::~ROAST_extendContigs_SCs() {
}

int max_coverage;
string ROAST_LE_tag = "trLE";

string ROAST_RE_tag = "trRE";

string ROAST_BE_tag = "trBE";

void ROAST_extendContigs_SCs::extend_bySoftclips(string fastafile, string ID_file, string BAM_FILE,  string num, int SC_support_TH, string tool_name_flag, int min_sc_reads, int sc_pos_from_corner, string exe_path) {

    utils utils;
    int base_quality = 0;
    std::size_t found = fastafile.find_last_of("/\\");
    string PATH = fastafile.substr(0, found);

    string bam_file = PATH + "/thread_" + num + ".bam";
   
    // string extract_bam_for_IDs = "cat " + ID_file + " | xargs samtools view -b " + BAM_FILE + " -o " + bam_file; // extract bam regions based on IDs generated earlier for each thread
    // string extract_bam_for_IDs = "intersectBed -abam " + BAM_FILE  + " -b " + ID_file + " > " + bam_file;

    string extract_bam_for_IDs = "samtools view " + BAM_FILE + " -M -L " + ID_file + " -b  -o " + bam_file; // + "  2> /dev/null";

    std::system(extract_bam_for_IDs.c_str());

    string cov_file = PATH + "/thread_" + num + ".cov"; // write for each fasta file
    string command_bam_to_cov = "samtools depth -aa " + bam_file + " > " + cov_file; // + "  2> /dev/null";
    std::system(command_bam_to_cov.c_str());

    map<string, map<int, int> > coverage = utils.coverage_data(cov_file);
    map<string, list<int> > repeat_contigs;
    /*do alignment first and calculate coverage*/

    string extendfasta = PATH + "/ext_" + num + ".fasta";
    string iteration_log = PATH + "/iteration_log_" + num + ".txt";
    string SC_id_file = PATH + "/SC_id_" + num + ".txt";
   

    ofstream extended_fasta;
    extended_fasta.open(extendfasta.c_str());

    ofstream log;
    log.open(iteration_log.c_str());

    ofstream SC_id;
    SC_id.open(SC_id_file.c_str());

    bool left = false, right = false, toleft = true;;

    map<string, list<int> >::iterator it_RC;
    list<int> repeat_con;

    string left_cons_seq = " ", right_cons_seq = " ", header, contig_name, fasta_seq;

    int extended_count = 0, deletion = 0, insertion = 0, count = 0 ;
    BamReader reader;

    vector<int> clipsize, read_pos, gen_pos;
    int prv_ref_id = 0, cur_ref_id, LSC_extended_size, RSC_extended_size, len_read_contig;
    bool pad, new_contig = false;
    
    vector <string> left_SC_seq, right_SC_seq, SC_seq, softclip_list;    
    vector <int>left_SC_base_pos, right_SC_base_pos;
    std::map<string, map<int, int> >::iterator itr_bowtie;

        
    BamTools::BamAlignment al;
    vector<CigarOp> cigar;
    BamTools::RefData ref_contig;


    // Opens a bam file
    if (!reader.Open(bam_file)) {
        cerr << "Could not open BAM file." << endl;
        exit(0);
    }

    // returns reference sequence entries
    vector<RefData> ref;
    ref = reader.GetReferenceData();
    // cout << "got bam data" << endl;

    while (reader.GetNextAlignment(al)) { //contig_name[ , , , ]
     
        if (al.RefID != -1 && !(al.Name.empty())) { // to avoid unaligned reads which gives -1 RefID


            if (count == 0) { // to get first entry for comparison to next entry
                prv_ref_id = al.RefID;
                count++;
            }
            ref_contig = ref.at(al.RefID);
            cur_ref_id = al.RefID;

            if (prv_ref_id == cur_ref_id) { // deal one contig at a time

                contig_name = ref_contig.RefName; // alignment trim sequence id after space so no check to remove -> [E, M, S, U]
               // cout << "contig name" << contig_name << endl;
                cigar = al.CigarData;

                prv_ref_id = cur_ref_id;

                if (new_contig == false) { // to consider only one time for cov extraction

                    itr_bowtie = coverage.find(contig_name);

                    unsigned currentMax = 0;

                    if (itr_bowtie != coverage.end()) {
                        map<int, int> &innermap = itr_bowtie->second;

                        for (map<int, int>::iterator it3 = innermap.begin(); it3 != innermap.end(); it3++) {
                            
                            if (it3->second > currentMax) {
                                currentMax = it3->second;
                            }
                        }
                        max_coverage = currentMax;
                    }
                    new_contig = true;
                    softclip_list.reserve(max_coverage);
                    left_SC_seq.reserve(max_coverage);
                    right_SC_seq.reserve(max_coverage);
                    left_SC_base_pos.reserve(max_coverage);
                    right_SC_base_pos.reserve(max_coverage);
                }
                  //get softclip data using bamtools
                if (al.GetSoftClips(clipsize, read_pos, gen_pos, pad) == true) {

                    len_read_contig = ref_contig.RefLength;
                    for (vector<CigarOp>::iterator itt = cigar.begin(); itt < cigar.end(); itt++) {
                        char type = itt->Type;
                        if (type == 'D')
                            deletion = deletion + itt->Length; //extract from read position in case deletion found
                    }
                    if (clipsize.size() == 1) { // if one sided clip , go for further checking

                        if ((gen_pos[0] <= sc_pos_from_corner) && (clipsize[0] == read_pos[0])){// && clipsize[0] >= min_sc_len)) {
                            int sc_start = 0, sc_end = 0;
                            string base_qualities = al.Qualities;
                            string sc_seq = "";

                            for (int sc_pos = 0; sc_pos <= read_pos[0]; sc_pos++) {

                                 int quality = int(base_qualities.at(sc_pos)) - 33; // ascii code - 33 (for phred score)

                                if (quality >= base_quality) { // base quality is not in the list of low base quality symbols
                                    sc_start = 0;
                                    sc_end = sc_pos;
                                    sc_seq = al.QueryBases.substr(0, sc_end);
                                } else
                                    break;
                            }
                            if (!(sc_seq.empty())) {
                                left_SC_seq.push_back(sc_seq);
                                left_SC_base_pos.push_back(gen_pos[0]);
                            }
                            sc_seq = "";

                            //left_SC_base_pos.push_back(gen_pos[0]);
                            // left_SC_seq.push_back(al.QueryBases.substr(0, read_pos[0]));

                        } else if ((gen_pos[0] >= len_read_contig - sc_pos_from_corner) && (clipsize[0] != read_pos[0])) { // && clipsize[0] >= min_sc_len)) {

                            int sc_start = 0, sc_end = 0;
                            string base_qualities = al.Qualities;
                            string sc_seq = "";
                            sc_start = read_pos[0] - deletion;

                            for (int sc_pos = read_pos[0] - deletion; sc_pos <= al.QueryBases.length() - 1; sc_pos++) {

                                int quality = int(base_qualities.at(sc_pos)) - 33; // ascii code - 33 (for phred score)

                                if (quality >= base_quality) { // base quality is not in the list of low base quality symbols
                                    sc_end = sc_pos;
                                    sc_seq = al.QueryBases.substr(sc_start, sc_end);
                                } else
                                    break;
                            }
                            if (!(sc_seq.empty())) {

                                right_SC_base_pos.push_back(gen_pos[0] + 1);
                                right_SC_seq.push_back(sc_seq);
                            }
                            sc_seq = "";

                            //  right_SC_base_pos.push_back(gen_pos[0] + 1);
                            // right_SC_seq.push_back(al.QueryBases.substr(read_pos[0] - deletion, al.QueryBases.length() - 1));
                        }

                    } else if (clipsize.size() > 1) // if two sided clips, check if opposite side clip is < = 3 
                    {
                        if ((gen_pos[0] <= sc_pos_from_corner) && (clipsize[0] == read_pos[0]) && clipsize[1] <= 3) { //left SC && clipsize[0] >= min_sc_len

                            int sc_start = 0, sc_end = 0;
                            string base_qualities = al.Qualities;
                            string sc_seq = "";

                            for (int sc_pos = 0; sc_pos <= read_pos[0]; sc_pos++) {

                                int quality = int(base_qualities.at(sc_pos)) - 33; // ascii code - 33 (for phred score)

                                if (quality >= base_quality) { // base quality is not in the list of low base quality symbols
                                    sc_start = 0;
                                    sc_end = sc_pos;
                                    sc_seq = al.QueryBases.substr(0, sc_end);
                                } else
                                    break;
                            }
                            if (!(sc_seq.empty())) {
                                left_SC_seq.push_back(sc_seq);
                                left_SC_base_pos.push_back(gen_pos[0]);
                            }
                            sc_seq = "";

                            //left_SC_base_pos.push_back(gen_pos[0]);
                            // left_SC_seq.push_back(al.QueryBases.substr(0, read_pos[0]));

                        } else if ((gen_pos[1] >= len_read_contig - sc_pos_from_corner) && (clipsize[1] != read_pos[1]) && clipsize[0] <= 3) { //right SC && clipsize[1] >= min_sc_len

                            int sc_start = 0, sc_end = 0;
                            string base_qualities = al.Qualities;
                            string sc_seq = "";
                            sc_start = read_pos[1] - deletion;

                            for (int sc_pos = read_pos[1] - deletion; sc_pos <= al.QueryBases.length() - 1; sc_pos++) {

                                 int quality = int(base_qualities.at(sc_pos)) - 33; // ascii code - 33 (for phred score)

                                if (quality >= base_quality) { // base quality is not in the list of low base quality symbols
                                    sc_end = sc_pos;
                                    sc_seq = al.QueryBases.substr(sc_start, sc_end);
                                } else
                                    break;
                            }
                            if (!(sc_seq.empty())) {

                                right_SC_base_pos.push_back(gen_pos[1] + 1);
                                right_SC_seq.push_back(sc_seq);
                            }
                            sc_seq = "";

                          //  right_SC_base_pos.push_back(gen_pos[1] + 1);
                          //  right_SC_seq.push_back(al.QueryBases.substr(read_pos[1] - deletion, al.QueryBases.length() - 1));

                        }
                    }
                    clipsize.clear();
                    read_pos.clear();
                    gen_pos.clear();
                    deletion = 0, insertion = 0;

                }

            } else {
               // cout << "got sc data for one contig" << endl;
                stringstream ss;
                string header = contig_name;
                
                 fasta_seq = utils.extract_fasta(contig_name, fastafile);
                it_RC = repeat_contigs.find(contig_name);

                if (!(left_SC_base_pos.empty()) && !(right_SC_base_pos.empty())) { //when both sides found

                    if (left_SC_base_pos.size() >= min_sc_reads) {
                        
                        softclip_list = extract_maxSC_block(left_SC_base_pos, left_SC_seq);
                        
                        int pos = SC_start_pos(left_SC_base_pos) + 1; //+1 because to extract coverage from itr_bowtie indexing start from 1 not 0 but in remaining indexing start from zero
                        /*extraction complete*/
                        if (softclip_list.size() >= min_sc_reads) {

                            map<int, int>::iterator itr_bowtie_left = itr_bowtie->second.find(pos);
                            if (itr_bowtie_left->second == 0) {
                                itr_bowtie_left->second = softclip_list.size(); // itr_bowtie->second.find(pos + 1); // because
                            }
                            //  cout << softclip_list.size() << " " << (itr_bowtie_left->second) << endl;

                            if (itr_bowtie_left->second == 0) // to through exception incase of abnormal coverage in coverage file
                                goto out_of_leftSC;

                            if ((float(softclip_list.size()) / float(itr_bowtie_left->second))*100 >= SC_support_TH) {
                                left_cons_seq = consensus_seq_CAP3(softclip_list, toleft, num, exe_path, PATH);//"left");

                                LSC_extended_size = left_cons_seq.size();

                                if (LSC_extended_size > 1) {//= min_sc_len) { //sc length should be greater than  //commented 25dec2020

                                    string fasta_temp = left_cons_seq + fasta_seq.substr(pos - 1, fasta_seq.size());
                                    
                                    if (it_RC == repeat_contigs.end()) { //contig is not already written
                                        repeat_con.push_back(fasta_temp.size());
                                        repeat_contigs[contig_name] = repeat_con;
                                    } else {
                                        repeat_con = it_RC->second;

                                        if (repeat_con.size() > 2) { // contig repeated more than three times 

                                            list<int>::reverse_iterator rit = repeat_con.rbegin();
                                            rit++; // to extract the second last value of list

                                            if (labs(*rit - fasta_temp.size()) > 5) { //difference in extension of current and previous - 1 is > 5 consider it 

                                                it_RC->second.push_back(fasta_temp.size());

                                            } else { // contig repeated more than 3 times, with no appearent change in size 

                                                left_cons_seq = "";
                                                goto out_of_leftSC;

                                            }
                                        } else 

                                            it_RC->second.push_back(fasta_temp.size());
                                    }
                                    fasta_seq = left_cons_seq + fasta_seq.substr(pos - 1, fasta_seq.size());
                                    ss << header << ";leftSC:" << LSC_extended_size;
                                    header = ss.str();
                                    repeat_con.clear();
                                    ss.str(string());
                                    left = true;
                                }
                            }
                        }
                    }
                    // else continue;
out_of_leftSC:
                    ;

                    softclip_list.clear();

                    if (right_SC_base_pos.size() >= min_sc_reads) {
                        softclip_list = extract_maxSC_block(right_SC_base_pos, right_SC_seq);
                        int pos = SC_start_pos(right_SC_base_pos) + 1;

                        if (softclip_list.size() >= min_sc_reads) {
                            map<int, int>::iterator itr_bowtie_right = itr_bowtie->second.find(pos);
                            // map<int, int>::iterator itr_bowtie_right = coverage_bowtie.find(pos);
                            
                            if (itr_bowtie_right->second == 0)
                                itr_bowtie_right->second = softclip_list.size(); //itr_bowtie->second.find(pos - 1);

                            //    cout << softclip_list.size() << " " << (itr_bowtie_right->second) << endl;
                            if (itr_bowtie_right->second == 0)
                                goto out_of_rightSC;

                            if (float((softclip_list.size()) / float(itr_bowtie_right->second))*100 >= SC_support_TH) {

                                right_cons_seq = consensus_seq_CAP3(softclip_list, !toleft , num, exe_path, PATH);//"right");
                                RSC_extended_size = right_cons_seq.size();   
                                
                                if (RSC_extended_size > 1) { // min_sc_len) { //commented 25dec2020
                                    string fasta_temp;

                                    if (left) 
                                        fasta_temp = fasta_seq.substr(0, pos - 2 + LSC_extended_size + 1) + right_cons_seq; //
                                    else
                                        fasta_temp = fasta_seq.substr(0, pos - 1) + right_cons_seq;

                                    if (it_RC == repeat_contigs.end()) { // if contig is not already written

                                        repeat_con.push_back(fasta_temp.size());
                                        repeat_contigs[contig_name] = repeat_con;

                                    } else {
                                        repeat_con = it_RC->second;
                                        
                                        if (repeat_con.size() > 2) {
                                            list<int>::reverse_iterator rit = repeat_con.rbegin();
                                            rit++; // to extract the second last value of list
                                            
                                            if (labs(*(rit) - fasta_temp.size()) > 5) { //difference in extension of current and previous - 1 is > 5 consider it 

                                                //cout << contig_name << " repeat" << endl;
                                                it_RC->second.push_back(fasta_temp.size());

                                            } else { // contig repeated more than 3 times, ignore

                                                right_cons_seq = "";
                                                goto out_of_rightSC;

                                            }
                                        } else 

                                            it_RC->second.push_back(fasta_temp.size());
                                    }
                                    repeat_con.clear();
                                    if (left)
                                        fasta_seq = fasta_seq.substr(0, pos - 2 + LSC_extended_size + 1) + right_cons_seq;
                                    else
                                        fasta_seq = fasta_seq.substr(0, pos - 1) + right_cons_seq;
                                    
                                    ss << header << ";rightSC:" << RSC_extended_size;
                                    header = ss.str();
                                    ss.str(string());
                                    log << header << endl;
                                    right = true;
                               }
                            }
                        }
                    }

out_of_rightSC:
                    ;
                    softclip_list.clear();
                    if (!right && left)
                        log << header << endl; // enter for left only

                } else if (left_SC_base_pos.size() >= min_sc_reads) {

                    softclip_list = extract_maxSC_block(left_SC_base_pos, left_SC_seq);
                    int pos = SC_start_pos(left_SC_base_pos) + 1;

                    if (softclip_list.size() >= min_sc_reads) {
                        // map<int, int>::iterator itr_bowtie_left = coverage_bowtie.find(pos);
                        map<int, int>::iterator itr_bowtie_left = itr_bowtie->second.find(pos);
                        if (itr_bowtie_left->second == 0) //if all scs at given position
                            itr_bowtie_left->second = softclip_list.size(); // itr_bowtie->second.find(pos + 1);

                        //cout << softclip_list.size() << " " << (itr_bowtie_left->second) << endl;
                        if (itr_bowtie_left->second == 0) // to through exception incase of abnormal coverage in coverage file
                            goto out_of_leftSC2;

                        if (float((softclip_list.size()) / float(itr_bowtie_left->second))*100 >= SC_support_TH) {

                            left_cons_seq = consensus_seq_CAP3(softclip_list, toleft, num, exe_path, PATH);//"left");
                            LSC_extended_size = left_cons_seq.size();

                           if ( LSC_extended_size > 1 ) {//= min_sc_len) { // commented 25dec2020

                                string fasta_temp = left_cons_seq + fasta_seq.substr(pos - 1, fasta_seq.size());

                                if (it_RC == repeat_contigs.end()) { //contig is not already written

                                    repeat_con.push_back(fasta_temp.size());
                                    repeat_contigs[contig_name] = repeat_con;

                                } else {

                                    repeat_con = it_RC->second;

                                    if (repeat_con.size() > 2) {

                                        list<int>::reverse_iterator rit = repeat_con.rbegin();
                                        rit++; // to extract the second last value of list

                                        if (labs(*rit - fasta_temp.size()) > 5) { //difference in extension of current and previous - 1 is > 5 consider it 

                                            //cout << contig_name << " repeat" << endl;
                                            it_RC->second.push_back(fasta_temp.size());

                                        } else { // contig repeated more than 3 times, ignore

                                            left_cons_seq = "";
                                            goto out_of_leftSC2;
                                        }
                                    } else 
                                        it_RC->second.push_back(fasta_temp.size());
                                }
                                repeat_con.clear();
                                fasta_seq = left_cons_seq + fasta_seq.substr(pos - 1, fasta_seq.size());
                                ss << header << ";leftSC:" << LSC_extended_size;
                                header = ss.str();
                                ss.str(string());
                                log << header << endl;
                                left = true;
                                //extended_count++;
                            }
                        }
                    }
                    //       else continue;
out_of_leftSC2:
                    ;
                    softclip_list.clear();

                } else if (right_SC_base_pos.size() >= min_sc_reads) { //(!(right_SC_base_pos.empty())) {

                    softclip_list = extract_maxSC_block(right_SC_base_pos, right_SC_seq);
                    int pos = SC_start_pos(right_SC_base_pos) + 1;

                    if (softclip_list.size() >= min_sc_reads) {
                        // map<int, int>::iterator itr_bowtie_right = coverage_bowtie.find(pos);
                        map<int, int>::iterator itr_bowtie_right = itr_bowtie->second.find(pos);
                        if (itr_bowtie_right->second == 0)

                            itr_bowtie_right->second = softclip_list.size(); //itr_bowtie->second.find(pos - 1);

                        if (itr_bowtie_right->second == 0) // to through exception incase of abnormal coverage in coverage file

                            goto out_of_rightSC2;

                        if ((float(softclip_list.size()) / float(itr_bowtie_right->second))*100 >= SC_support_TH) {

                            right_cons_seq = consensus_seq_CAP3(softclip_list, !toleft, num, exe_path, PATH);//"right");
                            RSC_extended_size = right_cons_seq.size();

                            if (RSC_extended_size > 1) { //= min_sc_len) { //commented 25dec2020

                                string fasta_temp = fasta_seq.substr(0, pos - 1) + right_cons_seq;

                                if (it_RC == repeat_contigs.end()) { // if contig is not already written

                                    repeat_con.push_back(fasta_temp.size());
                                    repeat_contigs[contig_name] = repeat_con;
                                } else {

                                    repeat_con = it_RC->second;

                                    if (repeat_con.size() > 2) {
                                        list<int>::reverse_iterator rit = repeat_con.rbegin();
                                        rit++; // to extract the second last value of list
                                        if (labs(*rit - fasta_temp.size()) > 5) { //difference in extension of current and previous - 1 is > 5 consider it 

                                            it_RC->second.push_back(fasta_temp.size());

                                        } else { // contig repeated more than 3 times, ignore

                                            right_cons_seq = "";
                                            goto out_of_rightSC2;
                                        }
                                    } else {

                                        it_RC->second.push_back(fasta_temp.size());
                                    }
                                }
                                repeat_con.clear();
                                fasta_seq = fasta_seq.substr(0, pos - 1) + right_cons_seq;
                                ss << header << ";rightSC:" << RSC_extended_size;
                                header = ss.str();
                                log << header << endl;
                                ss.str(string());
                                right = true;
                                //extended_count++;
                            }
                        }//
                    }
                    //        else continue;
out_of_rightSC2:
                    ;
                    softclip_list.clear();
                }

                size_t update_ind = contig_name.find(tool_name_flag);
                size_t update_ind2 = contig_name.rfind(tool_name_flag);

                if (left || right) { //if extended

                    if (update_ind != string::npos) //[] found
                    {
                        //string update_str = contig_name.substr(update_ind + tool_name_flag.length() + 1, update_ind2); // get everything within :ROAST_ _ROAST,
                        //string update_str = contig_name.substr(update_ind + tool_name_flag.length() + 1, update_ind2 - (update_ind + tool_name_flag.length()));

                        size_t find_BE = contig_name.find(ROAST_BE_tag);
                        size_t find_LE = contig_name.find(ROAST_LE_tag);
                        size_t find_RE = contig_name.find(ROAST_RE_tag);

                        if (left && right) {

                            if (find_BE != string::npos) {
                                //ignore
                            } else if (find_RE != string::npos) {

                                contig_name = contig_name.replace(find_RE, ROAST_BE_tag.size(), ROAST_BE_tag);

                            } else if (find_LE != string::npos) {

                                contig_name = contig_name.replace(find_LE, ROAST_BE_tag.size(), ROAST_BE_tag);

                            } else {

                                contig_name = contig_name.substr(0, contig_name.size() - tool_name_flag.length() - 1) + "-" + ROAST_BE_tag + "_" + tool_name_flag;
                                extended_count++; //not already extended -> to get unique count of each extended contig

                            }

                            SC_id << contig_name << "\t" << left_cons_seq << "\t" << right_cons_seq << endl;

                        } else if (left) {

                            if (find_BE != string::npos) {
                                //ignore
                            } else if (find_RE != string::npos) {

                                contig_name = contig_name.replace(find_RE, ROAST_BE_tag.size(), ROAST_BE_tag);

                            } else if (find_LE != string::npos) {
                                //ignore
                            } else {

                                contig_name = contig_name.substr(0, contig_name.size() - tool_name_flag.length() - 1) + "-" + ROAST_LE_tag + "_" + tool_name_flag;
                                extended_count++; //not already extended -> to get unique count of each extended contig
                            }

                            SC_id << contig_name << "\t" << "L" << "\t" << left_cons_seq << endl;

                        } else{// if (right == true) {

                            if (find_BE != string::npos) {
                                //ignore
                            } else if (find_RE != string::npos) {
                                //ignore
                            } else if (find_LE != string::npos) {

                                contig_name = contig_name.replace(find_LE, ROAST_BE_tag.size(), ROAST_BE_tag);

                            } else {

                                contig_name = contig_name.substr(0, contig_name.size() - tool_name_flag.length() - 1) + "-" + ROAST_RE_tag + "_" + tool_name_flag;
                                extended_count++; //not already extended -> to get unique count of each extended contig

                            }
                            SC_id << contig_name << "\t" << "R" << "\t" << right_cons_seq << endl;
                        }
                        //}
                    } else { // NO improvement yet                 

                        if (left && right) {

                            contig_name = contig_name + "_" + tool_name_flag + "_" + ROAST_BE_tag + "_" + tool_name_flag;
                            SC_id << contig_name << "\t" << left_cons_seq << "\t" << right_cons_seq << endl;

                        } else if (left) {

                            contig_name = contig_name + "_" + tool_name_flag + "_" + ROAST_LE_tag + "_" + tool_name_flag;
                            SC_id << contig_name << "\t" << "L" << "\t" << left_cons_seq << endl;

                        } else { // if (right == true) {

                            contig_name = contig_name + "_" + tool_name_flag + "_" + ROAST_RE_tag + "_" + tool_name_flag;
                            SC_id << contig_name << "\t" << "R" << "\t" << right_cons_seq << endl;
                        }
                        extended_count++;
                    }
                }
                extended_fasta << ">" << contig_name << endl << fasta_seq << endl; // write for any contig_name update above

                fasta_seq = "";

                //SC_id.insert(left_right(contig_name, left_cons_seq)); //insert a bunch of values
                //SC_id.insert(left_right(contig_name, right_cons_seq)); 
                
                left_SC_seq.clear();
                right_SC_seq.clear();
                SC_seq.clear();
                left_SC_base_pos.clear();
                right_SC_base_pos.clear();
                prv_ref_id = cur_ref_id;
                ss.str(string());
                header = "";
                left_cons_seq.clear();
                right_cons_seq.clear();
                left = false, right = false;
                contig_name = ref_contig.RefName;
                new_contig = false;
            }
        }
    }
     fasta_seq = utils.extract_fasta(contig_name, fastafile);
     extended_fasta << ">" << contig_name << endl << fasta_seq << endl; 
    // cout << "Extended incomplete contigs using soft-clipping at iteration " << super_iteration << "-" << mini_iteration << ": " << extended_count << endl;
    reader.Close();
    log.close();
    extended_fasta.close();
    SC_id.close();

    string thread_done = fastafile.substr(0, found)  + "/done_" + num + ".txt";
    ofstream done;
    done.open(thread_done.c_str());
    done << extended_count << endl;
    done.close();

}
vector<string> ROAST_extendContigs_SCs::extract_maxSC_block(vector<int> SC_base_pos, vector<string> SC_seq) {
    utils utils;
    int count = 1, max = 1, ind = 0, st_ind = 0, en_ind = 0, temp_ind = 0, pos = 0, pos_next;
    vector<string> softclip_list;
    softclip_list.reserve(max_coverage);
    vector<int> all_ind;
    all_ind.reserve(500);
    
    int prev = SC_base_pos[0]; 
    
    for (std::vector<int>::iterator it = SC_base_pos.begin(); it != SC_base_pos.end(); ++it) {
        
        if (prev == *it) {
            count++;
            pos_next = prev;

        } else {
            if (prev == pos) { // SC pos repeats 
                st_ind = temp_ind;
                temp_ind = ind; // + 1;
                en_ind = ind;
                if (count > max) max = count;
                // if (pos_next == pos) {
                all_ind.push_back(st_ind);
                all_ind.push_back(en_ind);
                pos = prev;
            } else if (count > max) { //new SC position
                st_ind = temp_ind;
                temp_ind = ind; // + 1;
                en_ind = ind;
                max = count;
                all_ind.clear();
                all_ind.push_back(st_ind);
                all_ind.push_back(en_ind);
                pos = prev;
                //  }
                // to get pattern like this {358 , 360, 360, 360, 359, 360, 360, 360 }

            } else
                temp_ind = ind; // + 1;
            count = 1;

        }
        ind++;
        prev = *it;
    }
    /*when only one position [0 0 0]*/
    if (count > max || prev == pos) {
        st_ind = temp_ind; //SC_base_pos.size() - count + 1;
        en_ind = SC_base_pos.size(); //st_ind + count -1 ;
        if (pos_next == pos) {
            all_ind.push_back(st_ind);
            all_ind.push_back(en_ind);
        } else if (count > max) {
            all_ind.clear();
            all_ind.push_back(st_ind);
            all_ind.push_back(en_ind);
        }
        pos = SC_base_pos[st_ind + 1];

    }
    if (!(all_ind.empty())) {
        for (int i = 0; i < all_ind.size() - 1; i = i + 2) { // to get the pattern like this {358 , 360, 360, 360, 359, 360, 360, 360 }
            st_ind = all_ind[i];
            en_ind = all_ind[i + 1];
            for (st_ind; st_ind < en_ind; st_ind++) {
                softclip_list.push_back(SC_seq[st_ind]);
            }
        }
    }

    return softclip_list; // soft-clipped reads block containing max SC at specific position       

}

int ROAST_extendContigs_SCs::SC_start_pos(vector<int> SC_base_pos) {
    utils utils;
    int count = 1, max = 1, ind = 0, st_ind = 0, en_ind = 0, temp_ind = 0, pos = 0;
    int prev = 0; // 
    
    for (std::vector<int>::iterator it = SC_base_pos.begin(); it != SC_base_pos.end(); ++it) {
        if (prev == *it) {
            count++;

        } else {
            if (count > max) {
                st_ind = temp_ind;
                temp_ind = ind + 1;
                en_ind = ind;
                max = count;
                pos = prev;
            } else
                temp_ind = ind + 1;
            count = 1;
        }
        ind++;
        prev = *it;
    }
    /*when only one position [0 0 0]*/
    if (count > max) {
        st_ind = SC_base_pos.size() - count;
        en_ind = st_ind + count;
        pos = SC_base_pos[st_ind + 1];
    }

    return pos; // soft-clipped reads block containing max SC at specific position         

}

struct length {

    bool operator()(const string& a, const string& b) {
        return a.size() < b.size();
    }
};

//void()

string ROAST_extendContigs_SCs::consensus_seq_CAP3(std::vector<string>softclip_list, bool toleft, string thread, string exe_path, string PATH) {

    string SCs_file = PATH + "/SCs_file_" + thread + ".fa";
    string CAP3_file = SCs_file + ".cap.contigs";

    ofstream SCs;
    SCs.open(SCs_file.c_str());
    
    int count = 1; vector<string> cons_seq_list;
    
    vector<string>::iterator it_SCs = softclip_list.begin(); // get SCs data from vector and make multi fasta file
    
    while(it_SCs != softclip_list.end()){
        
        SCs << ">" << count << endl << *it_SCs << endl;
        it_SCs++;
        count ++;
    }
    SCs.close();
    
    string cap3_assembly = exe_path + "external_tools/cap3 " + SCs_file + "  > /dev/null 2>&1";
    system(cap3_assembly.c_str());
    string temp = SCs_file + ".cap.ace";
    remove(temp.c_str());
    temp = SCs_file + ".cap.contigs.links";
    remove(temp.c_str());
    temp = SCs_file + ".cap.contigs.qual";
    remove(temp.c_str());
    temp = SCs_file + ".cap.info";
    remove(temp.c_str());
    temp = SCs_file + ".cap.singlets";
    remove(temp.c_str());

    ifstream assem;
    assem.open(CAP3_file.c_str());
    
    string cons_seq, entry;
    int cap3_assembly_count = 0;
    
    if (assem.peek() != std::ifstream::traits_type::eof()) { //empty file
        
        while (!(getline(assem, entry).eof())) { 

            if (entry[0] == '>') { //ignore header
                cons_seq_list.push_back(cons_seq);
                cons_seq = "";
                cap3_assembly_count++;
            } else
                cons_seq = cons_seq + entry;
        }
        assem.close();
        cons_seq_list.push_back(cons_seq);
        vector<string>::iterator it_cons_seq;
        it_cons_seq = cons_seq_list.begin();

        if (cap3_assembly_count == 1) { // for only one consensus sequences generated by CAP3 (ideal)
            it_cons_seq++;
            cons_seq = *it_cons_seq;
        } else { // get the longest consensus sequence from all generated by CAP3
           // cout << endl << "more than 1 consensus sequence found for extension using softclips " << thread << endl;
            sort(cons_seq_list.begin(), cons_seq_list.end(), length()); // sort list of SCs
            cons_seq = cons_seq_list.back(); // get consensus sequence of max length
        }
        //new_contig = utils.Rcomplement(new_contig); // take reverse complement of mates or new contig
    }
    if (!(cons_seq.empty())) {
        int cov_a = 0, consec_a = 0;
        int full_A_tail = 90;
        int partial_A_tail = 40;

        for (int i = 0; i < cons_seq.length(); i++) {
            if (cons_seq[i] == 'A') {
                cov_a++;
                consec_a++;
            } else
                consec_a = 0;
        }
        if ((float(cov_a / cons_seq.length())*100 >= full_A_tail) || (float(consec_a / cons_seq.length())*100 >= partial_A_tail)) {
            cons_seq = ""; // ignore extension
        }
    }
    else {  // if CAP3 doesn't generate any consensus sequennce (for shorter sequences), call consensus_seq function to generate consensus seq manually
         cons_seq = consensus_seq(softclip_list, toleft);
    }

remove(CAP3_file.c_str());
remove(SCs_file.c_str());
return cons_seq;

}


string ROAST_extendContigs_SCs::consensus_seq(std::vector<string>softclip_list, bool toleft) {
utils utils;
    float cov_a = 0, consec_a = 0; // for poly A tails
    float cov_c = 0;
    float cov_g = 0;
    float cov_t = 0;
    int size = 0;
    float count = 0;
    string final_cons_sequence;
    int count_amb_code = 0;
    int max_amb_code_allowed = 1;
    int full_A_tail = 90;
    int partial_A_tail = 40;
    int max_cov = 60; //75;
    string str = "";
    char code;
    sort(softclip_list.begin(), softclip_list.end(), length()); // sort list of SCs
    string s = softclip_list.back(); // get SCs of max length
    //string s = softclip_list[softclip_list.size()-1];
    size = s.size() - 1; // get size of max SC
    vector <char> con_seq; // [size-1] = "";
    con_seq.reserve(100);

    // cout << size << flag << endl;
    //  if (size > 1) {
    if (toleft) { //start

        int t = size;
        int clip_size = 1;
        while (t > 0) { // for each position till full length of max SC
            
            for (int z = 0; z < softclip_list.size(); z++) { // for each SC string
                
                str = softclip_list[z];
                int itt = str.size();
                if (itt <= clip_size)
                    continue;
                else {
                    if (str[itt - clip_size] == 'A')
                        cov_a = cov_a + 1;
                    else if (str[itt - clip_size] == 'C')
                        cov_c = cov_c + 1;
                    else if (str[itt - clip_size] == 'G')
                        cov_g = cov_g + 1;
                    else if (str[itt - clip_size] == 'T')
                        cov_t = cov_t + 1;
                    count++;
                }
                str = "";
            }

            if (count <= 2)
                break;
            else {
                if (float(cov_a / count)*100 >= max_cov) // (cov_a > cov_g && cov_a > cov_c && cov_a > cov_t)
                    con_seq.push_back('A');
                else if (float(cov_c / count)*100 >= max_cov) //(cov_c > cov_a && cov_c > cov_g && cov_c > cov_t)
                    con_seq.push_back('C');
                else if (float(cov_g / count)*100 >= max_cov) // (cov_g > cov_a && cov_g > cov_c && cov_g > cov_t)
                    con_seq.push_back('G');
                else if (float(cov_t / count)*100 >= max_cov) //(cov_t > cov_a && cov_t > cov_c && cov_t > cov_g)
                    con_seq.push_back('T');

                    // for two
                else if (count_amb_code < max_amb_code_allowed) {
                    
                    /*if (cov_c == cov_a && cov_c > cov_g && cov_c > cov_t) {
                        code = utils.base_code('C', 'A');
                        count_amb_code++;
                        con_seq.push_back(code);
                    } else if (cov_g == cov_a && cov_g > cov_c && cov_g > cov_t) {
                        code = utils.base_code('G', 'A');
                        count_amb_code++;
                        con_seq.push_back(code);
                    } else if (cov_t == cov_a && cov_t > cov_c && cov_t > cov_g) {
                        code = utils.base_code('T', 'A');
                        count_amb_code++;
                        con_seq.push_back(code);
                    } else if (cov_g == cov_c && cov_g > cov_a && cov_c > cov_t) {
                        code = utils.base_code('G', 'C');
                        count_amb_code++;
                        con_seq.push_back(code);
                    }//
                    else if (cov_g == cov_t && cov_g > cov_c && cov_g > cov_a) {
                        code = utils.base_code('G', 'T');
                        count_amb_code++;
                        con_seq.push_back(code);
                    } else if (cov_t == cov_c && cov_t > cov_a && cov_t > cov_g) {
                        code = utils.base_code('C', 'T');
                        count_amb_code++;
                        con_seq.push_back(code);
                    } else { */
                        code = 'N';
                        count_amb_code++;
                        con_seq.push_back(code);

                 //   }
                }//for three bases g
                    //                else if (cov_g == cov_c && cov_g == cov_a && cov_g > cov_t) {
                    //                    code = m_sc.base_code('G', 'C','A');
                    //                    count_amb_code++;
                    //                    con_seq.push_back(code);
                    //                }//
                    //                else if (cov_g == cov_t && cov_g == cov_a && cov_g > cov_c) {
                    //                    code = m_sc.base_code('G', 'T' ,'A');
                    //                    count_amb_code++;
                    //                    con_seq.push_back(code);
                    //                } else if (cov_g == cov_c && cov_g == cov_t && cov_g > cov_a) {
                    //                    code = m_sc.base_code('C', 'T', 'G');
                    //                    count_amb_code++;
                    //                    con_seq.push_back(code);
                    //                }
                    //                else if (cov_a == cov_c && cov_a == cov_t && cov_a > cov_g) {
                    //                    code = m_sc.base_code('C', 'T','A');
                    //                    count_amb_code++;
                    //                    con_seq.push_back(code);
                    //                }                   
                    //                else
                    //                { 
                    //                    con_seq.push_back('N');
                    //                    count_amb_code++;
                    //                }             
                else
                    break; //don't allow 1 base mismatches or N
                clip_size++;
                t--;
                cov_a = cov_c = cov_g = cov_t = 0;
                count = 0;
            }
        }
        // if (count_amb_code <= max_amb_code_allowed) {
        std::reverse(con_seq.begin(), con_seq.end());
        string cons_sequence(con_seq.begin(), con_seq.end());
        final_cons_sequence = cons_sequence;
        //  } else final_cons_sequence == "";
        con_seq.empty();
    } else { //flag: right soft clipping
        int clip_size = 0;
        while (clip_size < size) {
            for (int z = 0; z < softclip_list.size(); z++) {
                str = softclip_list[z];
                if (clip_size >= str.size())continue;
                else {
                    if (str[clip_size] == 'A')
                        cov_a = cov_a + 1;
                    else if (str[clip_size] == 'C')
                        cov_c = cov_c + 1;
                    else if (str[clip_size] == 'G')
                        cov_g = cov_g + 1;
                    else if (str[clip_size] == 'T')
                        cov_t = cov_t + 1;
                    count++;
                }
                str = "";
            }

            if (count <= 2)
                break;
            else {
                if (float(cov_a / count)*100 >= max_cov) // (cov_a > cov_g && cov_a > cov_c && cov_a > cov_t)
                    con_seq.push_back('A');
                else if (float(cov_c / count)*100 >= max_cov) //(cov_c > cov_a && cov_c > cov_g && cov_c > cov_t)
                    con_seq.push_back('C');
                else if (float(cov_g / count)*100 >= max_cov) // (cov_g > cov_a && cov_g > cov_c && cov_g > cov_t)
                    con_seq.push_back('G');
                else if (float(cov_t / count)*100 >= max_cov) //(cov_t > cov_a && cov_t > cov_c && cov_t > cov_g)
                    con_seq.push_back('T');

                    //for two bases
                else if (count_amb_code < max_amb_code_allowed) {
                  /*  if (cov_c == cov_a && cov_c > cov_g && cov_c > cov_t) {
                        code = utils.base_code('C', 'A');
                        count_amb_code++;
                        con_seq.push_back(code);
                    }//
                    else if (cov_g == cov_a && cov_g > cov_c && cov_g > cov_t) {
                        code = utils.base_code('G', 'A');
                        count_amb_code++;
                        con_seq.push_back(code);
                    } else if (cov_t == cov_a && cov_t > cov_c && cov_t > cov_g) {
                        code = utils.base_code('T', 'A');
                        count_amb_code++;
                        con_seq.push_back(code);
                    } else if (cov_g == cov_c && cov_g > cov_a && cov_c > cov_t) {
                        code = utils.base_code('G', 'C');
                        count_amb_code++;
                        con_seq.push_back(code);
                    }//
                    else if (cov_g == cov_t && cov_g > cov_c && cov_g > cov_a) {
                        code = utils.base_code('G', 'T');
                        count_amb_code++;
                        con_seq.push_back(code);
                    } else if (cov_t == cov_c && cov_t > cov_a && cov_t > cov_g) {
                        count_amb_code++;
                        code = utils.base_code('C', 'T');
                        con_seq.push_back(code);
                    } else {*/
                        count_amb_code++;
                        code = 'N';
                        con_seq.push_back(code);

                    //}
                }//for three bases g
                    //                else if (cov_g == cov_c == cov_a && cov_g > cov_t) {
                    //                    code = m_sc.base_code('G', 'C','A');
                    //                    count_amb_code++;
                    //                    con_seq.push_back(code);
                    //                }//
                    //                else if (cov_g == cov_t == cov_a && cov_g > cov_c) {
                    //                    code = m_sc.base_code('G', 'T' ,'A');
                    //                    count_amb_code++;
                    //                    con_seq.push_back(code);
                    //                } else if (cov_g == cov_c == cov_t && cov_g > cov_a) {
                    //                    code = m_sc.base_code('C', 'T', 'G');
                    //                    count_amb_code++;
                    //                    con_seq.push_back(code);
                    //                }
                    //                else if (cov_a == cov_c == cov_t && cov_a > cov_g) {
                    //                    code = m_sc.base_code('C', 'T','A');
                    //                    count_amb_code++;
                    //                    con_seq.push_back(code);
                    //                }              
                    //                else
                    //                { 
                    //                    con_seq.push_back('N');
                    //                    count_amb_code++;
                    //                }
                else
                    break; //don't allow 3 base mismatches or N
                clip_size++;
                cov_a = cov_c = cov_g = cov_t = 0;
                count = 0;
            }
        }
        //  if (count_amb_code <= max_amb_code_allowed) {
        string cons_sequence(con_seq.begin(), con_seq.end());
        final_cons_sequence = cons_sequence;
        //} else    final_cons_sequence = "";
        con_seq.empty();
        cov_a = 0;
        for (int i = 0; i < final_cons_sequence.length(); i++) {
            if (final_cons_sequence[i] == 'A') {
                cov_a++;
                consec_a++;
            } else consec_a = 0;
        }
        if ((float(cov_a / final_cons_sequence.length())*100 >= full_A_tail) || (float(consec_a / final_cons_sequence.length())*100 >= partial_A_tail))
            final_cons_sequence = ""; // ignore extension
    }

    //  }
    return final_cons_sequence;
}

int main(int argc, char** argv) {
    
 if (argc < 2) {
         cout << "No arguments found" << endl;
         exit(0);
     } 
      vector <string> args(argv, argv + argc);
      
      string fastafile = args[1];
      
      string ID_file = args[2];
      
      string bam_file = args[3];
            
      string num = args[4];
      
      int SC_support_TH = atoi(args[5].c_str()); //75
      
      string tool_name_flag = args[6];

      int min_sc_reads = atoi(args[7].c_str()); // 3
      
      int sc_pos_from_corner = atoi(args[8].c_str());
      
      string exe_path = args[9];
      
      ROAST_extendContigs_SCs tt;
      tt.extend_bySoftclips(fastafile, ID_file, bam_file, num, SC_support_TH, tool_name_flag, min_sc_reads, sc_pos_from_corner, exe_path);
      //cout <<"ROAST_extendContigs_SCs returned successfully " << endl;

      return 0;     
}
