/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   byRI.cpp
 * Author: madiha
 * 
 * Created on September 17, 2018, 1:53 PM
 */

#include "byRI.h"
#include "global.h"
#include "utils.h"
#include "ROAST_extendContigs.h"
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
#include <boost/regex.hpp>
#include <boost/filesystem/operations.hpp>

using namespace BamTools;
using namespace std;
using namespace boost;

byRI::byRI() {
}

byRI::byRI(const byRI& orig) {
}

byRI::~byRI() {
}

bool frgaments_merged = false;
int start_overlap, end_overlap;

struct RI_data {
    string read_contig;
    string mate_contig;
    int read_start;
    int read_end;
    bool strand_read; // reverse = true, forward = false
    int read_count;
    int mate_start; // start of first read
    int mate_end; // start of last read + read length  --> because we don't have end of mate 
    bool strand_mate; // reverse = true, forward = false //string strand_mate1;
    bool flag;

};

struct merged_info {
    string updated_id;
    string updated_seq;
    string updated_strand; // forward or RC for merging

};


int byRI::read_island(string bam_in, string in_file, string out_file, string &log, std::ofstream &log_time) {

    utils utils;
    RI_data ri_data;
    // create bam file containing only read islands and unmapped reads at the corners

    string bam_distant = path_inter + "/DISTANT.bam"; 
    string bam_unmapped = path_inter + "/UNMAPPED.bam";


    //ofstream temp_file;
    //temp_file.open(tempfile.c_str());
    utils.filter_sam_distant(bam_in, bam_distant);
    utils.filter_sam_unmapped(bam_in, bam_unmapped);

    string bam_indx_cmd = "samtools index " + bam_distant + " > /dev/null 2>&1";
    std::system(bam_indx_cmd.c_str());

    bam_indx_cmd = "samtools index " + bam_unmapped + " > /dev/null 2>&1";
    std::system(bam_indx_cmd.c_str());


    /*first extend assembly through CAP3 using unmapped reads*/
    time_t begin1, end1;
    time(&begin1);
    set<string> merged_check;


    ofstream log_cap3;
    log_cap3.open(log.c_str(), ios::app);
    log_cap3 << "outer iteration:" << outer_iteration << endl;
    log_cap3.close();

    //cout << "before cap3" << endl;
    extend_unmappedBam(bam_unmapped, in_file, out_file, log);
    //cout << "after cap3" << endl;
    /* Extension completed*/

    time(&end1);
    double elapsed_secs1 = double(end1 - begin1);
    log_time << "outer_iteration:" << outer_iteration << " CAP3 assembly process time:" << elapsed_secs1 << endl;

    ofstream mergedRI;
    mergedRI.open(out_file.c_str(), ios::app);
    ofstream log_new;
    log_new.open(log.c_str(), ios::app);

    int read_count = 0, count = 0, merged_count = 0;
    long prev_read_id = 0, prev_mate_id = 0;

    string entry, entry_m_sc, contig_RI, contig_MI, merged, RI_seq, MI_seq, header_fasta = "";

    float RI_start, MI_start, RI_end, MI_end, RI_IS, MI_IS, read_contig_len, mate_contig_len;

    multimap<int, RI_data> broken_byRI;
    typedef map<string, RI_data> temp_data;
    temp_data one_RI;

    string read_contig1, read_contig2, mate_contig1; // mate_contig2, strand_read1, strand_read2, strand_mate1, strand_mate2;
    bool strand_read1, strand_read2, strand_mate1, strand_mate2, RI_dir, MI_dir, flag; // reverse = true, forward = false

    int read_start1, read_end1, mate_start1, mate_end1, read_start2, read_end2, mate_start2, mate_end2;

    BamReader reader;
    const SamHeader header = reader.GetHeader(); // returns header data object
    BamTools::BamAlignment al;
    BamTools::RefData read_id, mate_id;
    vector<RefData> ref;

    if (!reader.Open(bam_distant)) {
        cerr << "Could not open BAM file." << endl;
        exit(0);
    }
    ref = reader.GetReferenceData();

    time_t start, end;
    time(&start);
    
    //cout << "before generating broken_byRI map" << endl;
    while (reader.GetNextAlignment(al)) {

        if (al.RefID != -1 && al.IsMateMapped() && al.IsMapped()) { // to avoid unaligned reads which gives -1 RefID

            read_id = ref.at(al.RefID);
            mate_id = ref.at(al.MateRefID);

            if (count == 0) {

                read_contig1 = read_id.RefName;
                mate_contig1 = mate_id.RefName;

                
                prev_read_id = al.RefID;
                prev_mate_id = al.MateRefID;

                read_start1 = al.Position;
                read_end1 = al.GetEndPosition();

                mate_start1 = al.MatePosition;
                mate_end1 = al.MatePosition  + read_length; //start of last read + read length
                
                strand_read1 = al.IsReverseStrand(); //reverse
                strand_mate1 = al.IsMateReverseStrand(); //reverse

                count += 1;

            } else {
                read_contig_len = read_id.RefLength;
                mate_contig_len = mate_id.RefLength;

                if (al.RefID == prev_read_id) { //if (read_contig1 == read_contig2) {

                    if (al.MateRefID == prev_mate_id) { //if (mate_contig1 == mate_contig2) {
                        read_end2 = al.GetEndPosition();

                       // if (read_count == 110)
                       //     cout << "wait" << endl;

                        if (al.MatePosition <= mate_start1) { // update mate start not mate end

                            mate_start2 = al.MatePosition; // to get start position of mate island
                            mate_end2 = mate_end1; // to get end position for mate island
                            mate_start1 = mate_start2; // to keep the first min value for start position

                        } else { // update mate end not mate start
                            mate_start2 = mate_start1;
                            mate_end2 = al.MatePosition + read_length;
                            mate_end1 = mate_end2; // to keep max end value for end position
                        }

                        strand_read2 = al.IsReverseStrand(); //true; // reverse
                        strand_mate2 = al.IsMateReverseStrand(); //true; //reverse



                        read_count++;

                    } else { // different mate contig
                       // read_count++;
                        temp_data::iterator it_temp = one_RI.find(mate_contig1);

                        if (it_temp != one_RI.end()) {

                            ri_data.read_count = it_temp->second.read_count + read_count;
                            one_RI.erase(it_temp);

                            if (read_start1 > it_temp->second.read_start)
                                read_start1 = it_temp->second.read_start;

                        } else {

                            ri_data.read_count = read_count;

                        }
                        ri_data.mate_start = mate_start1;
                        ri_data.read_start = read_start1;
                        ri_data.read_end = read_end1;
                        ri_data.mate_end = mate_end1;// + read_length; // to get end position of mate island
                        ri_data.strand_read = strand_read1;
                        ri_data.mate_contig = mate_contig1;
                        ri_data.read_contig = read_contig1;
                        ri_data.strand_mate = strand_mate1;
                        ri_data.flag = true; //active status

                        one_RI.insert(make_pair(mate_contig1, ri_data));
                        read_count = 0;
                        count = 0;
                    }
                } else {
                    if (read_count != 0) { // for last mate contig when read contig2 change instead mate_contig2  
                       // read_count++;
                        temp_data::iterator it_temp = one_RI.find(mate_contig1);

                        if (it_temp != one_RI.end()) {

                            ri_data.read_count = it_temp->second.read_count + read_count;
                            one_RI.erase(it_temp);

                        } else {

                            ri_data.read_count = read_count;

                        }
                        ri_data.mate_start = mate_start1;
                        ri_data.read_start = read_start1;
                        ri_data.read_end = read_end1;
                        ri_data.mate_end = mate_end1;// + read_length; // to get end position of mate island
                        ri_data.strand_read = strand_read1;
                        ri_data.mate_contig = mate_contig1;
                        ri_data.read_contig = read_contig1;
                        ri_data.strand_mate = strand_mate1;
                        ri_data.flag = true; //active status

                        one_RI.insert(make_pair(mate_contig1, ri_data));
                        read_count = 0;
                    }

                    for (temp_data::iterator itr = one_RI.begin(); itr != one_RI.end(); itr++) {
                        broken_byRI.insert(make_pair(itr->second.read_count, itr->second));
                        //temp_file << itr->second.read_count << " " << itr->second.read_contig << " " << itr->second.read_start << " " << itr->second.read_end << " " << itr->second.mate_contig << " " << itr->second.mate_start << " " << itr->second.mate_end << endl;
                    }
                    one_RI.clear();
                    read_count = 0;
                    count = 0;
                }
                read_contig1 = read_id.RefName;
                mate_contig1 = mate_id.RefName;
                prev_read_id = al.RefID;
                prev_mate_id = al.MateRefID;
                read_end1 = read_end2;
                strand_read1 = strand_read2;
                strand_mate1 = strand_mate2;
            }
        }
    }// end of while

    for (temp_data::iterator itr = one_RI.begin(); itr != one_RI.end(); itr++) {
        broken_byRI.insert(make_pair(itr->second.read_count, itr->second));
        //temp_file << itr->second.read_count << " " << itr->second.read_contig << " " << itr->second.read_start << " " << itr->second.read_end << " " << itr->second.mate_contig << " " << itr->second.mate_start << " " << itr->second.mate_end << endl;
    }
    one_RI.clear();
    read_count = 0;

    //cout << "after generating broken_byRI map" << endl;

    // broken_byRI is in sorted order per following rules
    // saved in ascending sorted order by read count
    // for same outer key it get the first one, otherwise it sort based on key of inner map
    // and if inner keys are same then it sort based on outer keys..

    typedef map<string, string> merged_fragment; // entries of this vectors are dir(read contig), updated id, updated seq, also separate for mate_contig
    merged_fragment merged_fragments;

    typedef map<string, vector <string> > id_mapp;
    id_mapp id_map;

    typedef map<string, string> final_merged_frag;
    final_merged_frag final_merged_frags;
    set<string> dir_update;

    multimap<int, RI_data>::reverse_iterator it_RI;
    multimap<int, RI_data>::iterator it_data;
    string new_seq;
    map<string, string>::iterator it_fasta;

    for (it_RI = broken_byRI.rbegin(); it_RI != broken_byRI.rend(); ++it_RI) { //map<sting, map <int, RI_data> >

        ri_data = it_RI->second;

        contig_RI = ri_data.read_contig;
        contig_MI = ri_data.mate_contig;
        RI_start = ri_data.read_start;
        MI_start = ri_data.mate_start;
        RI_end = ri_data.read_end;
        MI_end = ri_data.mate_end;
        RI_dir = ri_data.strand_read;
        MI_dir = ri_data.strand_mate;
        flag = ri_data.flag;
        read_count = ri_data.read_count;

        it_fasta = AllFasta_data.find(contig_RI);
        RI_seq = it_fasta->second;
        
        it_fasta = AllFasta_data.find(contig_MI);
        MI_seq = it_fasta->second;
                
       // RI_seq = utils.extract_fasta(contig_RI, in_file);
       // MI_seq = utils.extract_fasta(contig_MI, in_file);
        
       // cout << contig_RI << " and " << contig_MI << endl;
        
       // if(contig_RI == "TART_57820_tr_trBE_tr" || contig_MI == "TART_57820_tr_trBE_tr")
        //    cout << "debugg" << endl;
         
        if (extended_byCAP3.find(contig_RI) == extended_byCAP3.end() && extended_byCAP3.find(contig_MI) == extended_byCAP3.end() && read_count >= min_distMapped_reads) {
            string idd1 = contig_RI;
            string idd2 = contig_MI;

            size_t update_ind = contig_RI.find(tool_name_flag); //tr

            if (update_ind != std::string::npos) { // [] found remove it and update 
                idd1 = contig_RI.substr(0, update_ind - 1);
            }

            update_ind = contig_MI.find(tool_name_flag);

            if (update_ind != std::string::npos) { // [] found remove it
                idd2 = contig_MI.substr(0, update_ind - 1);
            }

            header_fasta = idd1;

            if (flag) {

                string ch1 = contig_RI;
                string ch2 = contig_MI;

                string ch3 = contig_MI;
                string ch4 = contig_RI;
                
                size_t chimera11 = contig_RI.rfind(ROAST_newconitg1); // because multi-gene Chimera has _a and _b at end of their ids before tag
                size_t chimera12 = contig_MI.rfind(ROAST_newconitg2);

                size_t chimera21 = contig_MI.rfind(ROAST_newconitg1); // because multi-gene Chimera has _a and _b at end of their ids before tag
                size_t chimera22 = contig_RI.rfind(ROAST_newconitg2);

                //   if (contig_RI != contig_MI) { have already check in the start
                if(chimera11  != string::npos)
                    ch1 = contig_RI.substr(0, chimera11 - 1);
                
                if(chimera12  != string::npos)
                    ch2 = contig_MI.substr(0, chimera12 - 1);
                
                if(chimera21  != string::npos)
                    ch3 = contig_MI.substr(0, chimera21 - 1);
                
                if(chimera22  != string::npos)
                    ch4 = contig_RI.substr(0, chimera22 - 1);

                if ((ch1 != ch2) && (ch4 != ch3)) { //either both contigs are not same or if they are chimeras split by _a, _b then their original contig ids should not same before _a and _b
                    // merge them
                    // find read contig if read strand is same change flag to false
                    // find mate contig (which appears as read now)if in same direction as just dealt here make flag false
                    // we also need read strand information and only change read contig direction to make it trackable


                    //check if contig_RI or contig_MI is already updated? use their updated name, seq, and strands

                    map<string, string> ::iterator it_ri;
                    it_ri = merged_fragments.find(contig_RI);

                    map<string, string> ::iterator it_mi;
                    it_mi = merged_fragments.find(contig_MI);

                    if (it_ri != merged_fragments.end()) {
                        //cout << "Read contig is repeating " << endl;
                        
                        if (dir_update.find(contig_RI) != dir_update.end()) {

                            if (!RI_dir) //forward
                            { 
                                int st = RI_start;
                                RI_start = RI_seq.size() - RI_end; // 9july21, to consider new direction after previous merging change start and end positions as well
                                RI_end = RI_seq.size() - st;
                                
                                RI_dir = true; //reverse
                            
                            }
                            else{ //reverse
                                int st = RI_start;
                                RI_dir = false;
                                RI_start = RI_seq.size() - RI_end; // 9july21, to consider new direction after previous merging change start and end positions as well
                                RI_end = RI_seq.size() - st;
                            }
                        }
                        //
                        header_fasta = it_ri->second; //

                        merged_count = merged_count + merge_RIs(final_merged_frags, id_map, merged_fragments, contig_RI, contig_MI, dir_update, RI_seq, MI_seq, RI_start, RI_end, MI_start, MI_end, RI_dir, MI_dir, header_fasta, idd2, it_ri, log_new, mergedRI, merged_check);

                    } else if (it_mi != merged_fragments.end()) {

                        //cout << "Mate contig is repeating" << endl;
                        if (dir_update.find(contig_MI) != dir_update.end()) { // if MI already done make it RI

                            if (!MI_dir) //forward
                            {
                                int st = MI_start;
                                MI_start = MI_seq.size() - MI_end; // 9july21, to consider new direction after previous merging change start and end positions as well
                                MI_end = MI_seq.size() - st;

                                MI_dir = true; //reverse

                            } else { //reverse

                                MI_dir = false;
                                int st = MI_start;
                                MI_start = MI_seq.size() - MI_end; // 9july21, to consider new direction after previous merging change start and end positions as well
                                MI_end = MI_seq.size() - st;
                            }
                        }

                        header_fasta = it_mi->second; //

                        merged_count = merged_count + merge_RIs(final_merged_frags, id_map, merged_fragments, contig_MI, contig_RI, dir_update, MI_seq, RI_seq, MI_start, MI_end, RI_start, RI_end, MI_dir, RI_dir, header_fasta, idd2, it_mi, log_new, mergedRI, merged_check);

                    } else { // both not already written

                       // cout << "new read and mate contigs" << endl;
                        if (!RI_dir) { //forward

                            if (MI_dir) { //reverse
                           //   cout << "FR" << endl;

                                merged = overlap_merge(RI_seq, MI_seq, RI_start, RI_end, MI_start, MI_end ); //, RI_pos);

                                if (!(merged.empty())) // no overlap found so add new batch number to both contigs which shows read island but can't be merged
                                {
                                    if ((contig_RI.find(ROASTcap3_left_tag) != std::string::npos) || (contig_MI.find(ROASTcap3_left_tag) != std::string::npos) || (contig_RI.find(ROASTcap3_right_tag) != std::string::npos) || (contig_MI.find(ROASTcap3_right_tag) != std::string::npos)) {

                                        int init_size; 
                                        bool RI_largest;

                                        if (RI_seq.size() > MI_seq.size()){
                                            init_size = RI_seq.size();
                                            RI_largest = true;
                                        }
                                        else{
                                            init_size = MI_seq.size();
                                            RI_largest = false;
                                        }
                                        if (merged.size() > init_size  && (merged.size() - init_size >= (min_CAP3_ext * read_length) / 100)) {
                                            
                                            header_fasta = header_fasta + "_and_" + idd2;

                                            merged_count++;

                                            // merged_contigs.push_back(header_fasta);
                                            // 28oct21 // write for chimera
                                            merged_contigs.insert(make_pair(header_fasta, start_overlap));
                                            merged_contigs.insert(make_pair(header_fasta, end_overlap));

                                            log_new << contig_RI << "\t" << contig_MI << "\t" << "RI(FR)," << "\t" << "overlap=" << start_overlap << "-" << end_overlap << endl;
                                        } else {
                                            merged_check.insert(contig_RI); // so that it won't get read again to write
                                            merged_check.insert(contig_MI);
                                            
                                            if (RI_largest) { // tag RI seq, discard MI

                                                string contig_RI_new = contig_RI + "_" + stop_RCAP3_ext;
                                                mergedRI << ">" << contig_RI_new << endl << RI_seq << endl;

                                            } else {

                                                string contig_MI_new = contig_MI + "_" + stop_LCAP3_ext;
                                                mergedRI << ">" << contig_MI_new << endl << MI_seq << endl;
                                            }
                                            merged.clear();

                                        }
                                    } else {
                                        string temp_header = header_fasta + "_and_" + idd2;
                                        
                                        if(temp_header.length() <= max_header_size ) //27dec 21  to avoid error in cufflink due to long ID size
                                            header_fasta = temp_header ;
                                        //else keep header_fasta as it is 

                                        merged_count++;

                                        // merged_contigs.push_back(header_fasta);
                                         // 28oct21 // write for chimera
                                            merged_contigs.insert(make_pair(header_fasta, start_overlap));
                                            merged_contigs.insert(make_pair(header_fasta, end_overlap));

                                        log_new << contig_RI << "\t" << contig_MI << "\t" << "RI(FR)," << "\t" << "overlap=" << start_overlap << "-" << end_overlap << endl;
                                    }
                                }


                            } else if (!MI_dir) { // forward/ take reverse complement
                                 // cout << "FF" << endl;
                                
                                string RI_seq_RC = utils.Rcomplement(RI_seq); //now its R

                                int ri_end = RI_seq.size() - RI_start;
                                int ri_start = RI_seq.size() - RI_end;

                                merged = overlap_merge(MI_seq, RI_seq_RC, MI_start, MI_end, ri_start, ri_end); //, MI_pos);
                                
                                if (!(merged.empty())) {
                                    
                                    if ((contig_RI.find(ROASTcap3_left_tag) != std::string::npos) || (contig_MI.find(ROASTcap3_left_tag) != std::string::npos) || (contig_RI.find(ROASTcap3_right_tag) != std::string::npos) || (contig_MI.find(ROASTcap3_right_tag) != std::string::npos)) {

                                        int init_size;
                                        bool RI_largest;

                                        if (RI_seq.size() > MI_seq.size()) {
                                            init_size = RI_seq.size();
                                            RI_largest = true;
                                        } else {
                                            init_size = MI_seq.size();
                                            RI_largest = false;
                                        }
                                        if (merged.size() > init_size  && ( merged.size() - init_size >= (min_CAP3_ext * read_length) / 100) ) {

                                           string temp_header = header_fasta + "_and_" + idd2;
                                        
                                        if(temp_header.length() <= max_header_size ) //27dec 21  to avoid error in cufflink due to long ID size
                                            header_fasta = temp_header ;
                                        //else keep header_fasta as it is 

                                            merged_count++;

                                            // merged_contigs.push_back(header_fasta);
                                             // 28oct21 // write for chimera
                                            merged_contigs.insert(make_pair(header_fasta, start_overlap));
                                            merged_contigs.insert(make_pair(header_fasta, end_overlap));  
                                            dir_update.insert(contig_RI);

                                            log_new << contig_RI << "\t" << contig_MI << "\t" << "RI(FF)," << "\t" << "overlap=" << start_overlap << "-" << end_overlap << endl;
                                        } else {
                                            merged_check.insert(contig_RI); // so that it won't get read again to write
                                            merged_check.insert(contig_MI);
                                            
                                          if (RI_largest) { // tag RI seq, discard MI

                                                string contig_RI_new = contig_RI + "_" + stop_RCAP3_ext;
                                                mergedRI << ">" << contig_RI_new << endl << RI_seq << endl;

                                            } else {

                                                string contig_MI_new = contig_MI + "_" + stop_LCAP3_ext;
                                                mergedRI << ">" << contig_MI_new << endl << MI_seq << endl;
                                            }
                                            merged.clear();

                                        }
                                    } else {
                                        string temp_header = header_fasta + "_and_" + idd2;
                                        
                                        if(temp_header.length() <= max_header_size ) //27dec 21  to avoid error in cufflink due to long ID size
                                            header_fasta = temp_header ;
                                        //else keep header_fasta as it is 

                                        merged_count++;

                                        // merged_contigs.push_back(header_fasta);
                                             // 28oct21 // write for chimera
                                            merged_contigs.insert(make_pair(header_fasta, start_overlap));
                                            merged_contigs.insert(make_pair(header_fasta, end_overlap)); 
                                        dir_update.insert(contig_RI);

                                        log_new << contig_RI << "\t" << contig_MI << "\t" << "RI(FF)," << "\t" << "overlap=" << start_overlap << "-" << end_overlap << endl;
                                    }
                                }

                            }


                        } else if (RI_dir) { //R

                            if (MI_dir) { //R
                              //	cout << "RR" << endl;
                                string RI_seq_RC = utils.Rcomplement(RI_seq); //now its F

                                int ri_end = RI_seq.size() - RI_start;
                                int ri_start = RI_seq.size() - RI_end;

                                merged = overlap_merge(RI_seq_RC, MI_seq, ri_start, ri_end, MI_start, MI_end); //, ri_pos);

                                if (!(merged.empty())) {

                                    if ((contig_RI.find(ROASTcap3_left_tag) != std::string::npos) || (contig_MI.find(ROASTcap3_left_tag) != std::string::npos) || (contig_RI.find(ROASTcap3_right_tag) != std::string::npos) || (contig_MI.find(ROASTcap3_right_tag) != std::string::npos)) {

                                        int init_size;
                                        bool RI_largest;

                                        if (RI_seq.size() > MI_seq.size()) {
                                            init_size = RI_seq.size();
                                            RI_largest = true;
                                        } else {
                                            init_size = MI_seq.size();
                                            RI_largest = false;
                                        }
                                        if (merged.size() > init_size && (merged.size() - init_size >= (min_CAP3_ext * read_length) / 100)) {
                                            string temp_header = header_fasta + "_and_" + idd2;

                                            if (temp_header.length() <= max_header_size) //27dec 21  to avoid error in cufflink due to long ID size
                                                header_fasta = temp_header;
                                            //else keep header_fasta as it is 

                                            merged_count++;

                                            // merged_contigs.push_back(header_fasta);
                                            // 28oct21 // write for chimera
                                            merged_contigs.insert(make_pair(header_fasta, start_overlap));
                                            merged_contigs.insert(make_pair(header_fasta, end_overlap));
                                            dir_update.insert(contig_RI);

                                            log_new << contig_RI << "\t" << contig_MI << "\t" << "RI(RR)," << "\t" << "overlap=" << start_overlap << "-" << end_overlap << endl;
                                            
                                        } else {
                                            merged_check.insert(contig_RI); // so that it won't get read again to write
                                            merged_check.insert(contig_MI);
                                         if (RI_largest) { // tag RI seq, discard MI

                                                string contig_RI_new = contig_RI + "_" + stop_RCAP3_ext;
                                                mergedRI << ">" << contig_RI_new << endl << RI_seq << endl;

                                            } else {

                                                string contig_MI_new = contig_MI + "_" + stop_LCAP3_ext;
                                                mergedRI << ">" << contig_MI_new << endl << MI_seq << endl;
                                            }
                                            merged.clear();

                                        }
                                    } else {

                                   string temp_header = header_fasta + "_and_" + idd2;
                                        
                                        if(temp_header.length() <= max_header_size ) //27dec 21  to avoid error in cufflink due to long ID size
                                            header_fasta = temp_header ;
                                        //else keep header_fasta as it is 

                                        merged_count++;

                                        // merged_contigs.push_back(header_fasta);
                                        // 28oct21 // write for chimera
                                        merged_contigs.insert(make_pair(header_fasta, start_overlap));
                                        merged_contigs.insert(make_pair(header_fasta, end_overlap));
                                        dir_update.insert(contig_RI);

                                        log_new << contig_RI << "\t" << contig_MI << "\t" << "RI(RR)," << "\t" << "overlap=" << start_overlap << "-" << end_overlap << endl;
                                    }
                                }

                            } else if (!MI_dir) { //F
                              //  cout << "RF" << endl;
                                merged = overlap_merge(MI_seq, RI_seq, MI_start, MI_end, RI_start, RI_end); //, MI_pos);

                                if (!(merged.empty())) {
                                    
                                    if ((contig_RI.find(ROASTcap3_left_tag) != std::string::npos) || (contig_MI.find(ROASTcap3_left_tag) != std::string::npos) || (contig_RI.find(ROASTcap3_right_tag) != std::string::npos) || (contig_MI.find(ROASTcap3_right_tag) != std::string::npos)) {

                                        int init_size;
                                        bool RI_largest;

                                        if (RI_seq.size() > MI_seq.size()) {
                                            init_size = RI_seq.size();
                                            RI_largest = true;
                                        } else {
                                            init_size = MI_seq.size();
                                            RI_largest = false;
                                        }
                                        if (merged.size() > init_size && (merged.size() - init_size >= (min_CAP3_ext * read_length) / 100)) {

                                         string temp_header = header_fasta + "_and_" + idd2;
                                        
                                        if(temp_header.length() <= max_header_size ) //27dec 21  to avoid error in cufflink due to long ID size
                                            header_fasta = temp_header ;
                                        //else keep header_fasta as it is 

                                            merged_count++;
                                            // merged_contigs.push_back(header_fasta);
                                            // 28oct21 // write for chimera
                                            merged_contigs.insert(make_pair(header_fasta, start_overlap));
                                            merged_contigs.insert(make_pair(header_fasta, end_overlap));
                                            log_new << contig_RI << "\t" << contig_MI << "\t" << "RI(RF)," << "\t" << "overlap=" << start_overlap << "-" << end_overlap << endl;

                                        } else {
                                            merged_check.insert(contig_RI); // so that it won't get read again to write 
                                            merged_check.insert(contig_MI);
                                           if (RI_largest) { // tag RI seq, discard MI

                                                string contig_RI_new = contig_RI + "_" + stop_RCAP3_ext;
                                                mergedRI << ">" << contig_RI_new << endl << RI_seq << endl;

                                            } else {

                                                string contig_MI_new = contig_MI + "_" + stop_LCAP3_ext;
                                                mergedRI << ">" << contig_MI_new << endl << MI_seq << endl;
                                            }
                                            merged.clear();
                                        }
                                    } else {
                                          string temp_header = header_fasta + "_and_" + idd2;
                                        
                                        if(temp_header.length() <= max_header_size ) //27dec 21  to avoid error in cufflink due to long ID size
                                            header_fasta = temp_header ;
                                        //else keep header_fasta as it is 
                                        merged_count++;
                                          // merged_contigs.push_back(header_fasta);
                                             // 28oct21 // write for chimera
                                            merged_contigs.insert(make_pair(header_fasta, start_overlap));
                                            merged_contigs.insert(make_pair(header_fasta, end_overlap)); 
                                        log_new << contig_RI << "\t" << contig_MI << "\t" << "RI(RF)," << "\t" << "overlap=" << start_overlap << "-" << end_overlap << endl;
                                    }
                                }
                            }
                        }
                        if (!(merged.empty())) { //fill all four data structures
                            merged_fragments.insert(make_pair(contig_RI, header_fasta));
                            merged_fragments.insert(make_pair(contig_MI, header_fasta));

                            vector <string> temp; // add new and vector of old ids in id map
                            temp.push_back(contig_RI);
                            temp.push_back(contig_MI);
                            id_map.insert(make_pair(header_fasta, temp));
                            temp.clear();

                            final_merged_frags.insert(make_pair(header_fasta, merged)); // add new id and new sequence
                            merged_check.insert(contig_RI);
                            merged_check.insert(contig_MI);
                            frgaments_merged = true;
                            
                            if(merged.size() <= RI_seq.size()) // to check if size of contig didn't increase after merging of fragmented contigs - > embeded sequences
                                embed_count++ ;
                        }
                        merged.clear();

                    }

                    //considering broken_byRI is sorted by KEY which is read island count

                    // for current and all new entries if same contig RI appear with same strand? ignore it  
                    if (frgaments_merged) {
                        
                        for (it_data = broken_byRI.begin(); it_data != broken_byRI.end(); it_data++) {

                            if (it_data->second.read_contig == contig_RI && it_data->second.strand_read == RI_dir)
                                it_data->second.flag = false;
                            if (it_data->second.read_contig == contig_MI && it_data->second.strand_read == MI_dir)
                                it_data->second.flag = false;

                            if (it_data->second.mate_contig == contig_RI && it_data->second.strand_mate == RI_dir)
                                it_data->second.flag = false;
                            if (it_data->second.mate_contig == contig_MI && it_data->second.strand_mate == MI_dir)
                                it_data->second.flag = false;
                        }
                    }

                    frgaments_merged = false;
                    header_fasta = "";
                }
                //}

                // find make all flags false here containing same contig_RI and direction and contig_MI and direction in broken_byRI map<sting, map <int, RI_data> >
            }
        }// cap3 check
    }

    for (final_merged_frag::iterator it_merged = final_merged_frags.begin(); it_merged != final_merged_frags.end(); it_merged++) {
        mergedRI << ">" << it_merged->first << endl << it_merged->second << endl;
        //cout << ">" << it_merged->first << endl << it_merged->second << endl;
    }

    time(&end);
    double elapsed_secs = double(end - start);
    log_time << "outer_iteration:" << outer_iteration << " Read islands process time:" << elapsed_secs << endl << endl  << "Embeded fragmented contigs: " << embed_count << endl;


    ifstream infile;
    infile.open(in_file.c_str());
    string ID, ID1;

    while (!(getline(infile, entry_m_sc).eof())) { // rewrite merge file for RI and keep SC merged as it is
        
        if (entry_m_sc[0] == '>') {
            ID = entry_m_sc;
            ID1 = entry_m_sc.substr(1, entry_m_sc.size());

        } else {
            if ((merged_check.find(ID1) == merged_check.end()) && (extended_byCAP3.find(ID1) == extended_byCAP3.end())) {//current contig is not being merged
                mergedRI << ID << endl << entry_m_sc << endl;
            } else {
                continue;
            }
        }
    }

    cout << "Merged fragmented contigs using Reads Islands at iteration " << outer_iteration << ": " << merged_count << endl << endl;
    merged_fragmented_contigs = merged_fragmented_contigs + merged_count; // GLOBAL COUNT
    infile.close();
    mergedRI.close();
    log_new.close();
    extended_byCAP3.clear();
    //temp_file.close();

    return merged_count;
}

int byRI::merge_RIs(map<string, string> &final_merged_frags, map<string, vector <string> > &id_map, map<string, string> &merged_fragments, string contig_RI2, string contig_MI2,
        set<string> &dir_update, string RI_seq2, string MI_seq2, int RI_start2, int RI_end2, int MI_start2, int MI_end2, bool RI_dir2, bool MI_dir2, string header_fasta, string idd2, map<string, string> ::iterator &it_ri, ofstream &log_new, ofstream &mergedRI, set <string> &merged_check) 
{
    utils utils;
    typedef map<string, string> merged_fragment; // entries of this vectors are dir(read contig), updated id, updated seq, also separate for mate_contig
    typedef map<string, vector <string> > id_mapp;
    typedef map<string, string> final_merged_frag;

    multimap<int, RI_data>::reverse_iterator it_RI;
    string new_seq, merged, new_header_fasta;

    int merged_count = 0;

    final_merged_frag::iterator it_update_seq;
    it_update_seq = final_merged_frags.find(header_fasta);

    string old_seq = RI_seq2; // 
    RI_seq2 = it_update_seq->second; // use previously merged updated seq for found contig_RI            

  //   cout << old_seq << endl  <<"old seq size " << old_seq.size() << endl;
  //   cout << RI_seq2 << endl << "RI size " << RI_seq2.size() << endl;
   //  cout << MI_seq2 << endl << "MI size" << MI_seq2.size() << endl;

    //update RI_pos here 

    id_mapp::iterator it_newid;
    it_newid = id_map.find(header_fasta);

    if (!RI_dir2) { //F

        if (MI_dir2) { //R
         //    cout << "FR" << endl;
            // MI is on right means RI has something on its left earlier 213 
            RI_start2 = RI_seq2.size() - (old_seq.size() - RI_start2); //update position 
            RI_end2 = RI_seq2.size() - (old_seq.size() - RI_end2);

            merged = overlap_merge(RI_seq2, MI_seq2, RI_start2, RI_end2, MI_start2, MI_end2); //, RI_pos); // RI_seq - RI_pos 

            if (!(merged.empty())) {

                if ((contig_RI2.find(ROASTcap3_left_tag) != std::string::npos) || (contig_MI2.find(ROASTcap3_left_tag) != std::string::npos) || (contig_RI2.find(ROASTcap3_right_tag) != std::string::npos) || (contig_MI2.find(ROASTcap3_right_tag) != std::string::npos)) {

                    int init_size;
                    bool RI_largest;

                    if (RI_seq2.size() > MI_seq2.size()) {
                        init_size = RI_seq2.size();
                        RI_largest = true;
                    } else {
                        init_size = MI_seq2.size();
                        RI_largest = false;
                    }
                    if (merged.size() > init_size  && (merged.size() - init_size >= (min_CAP3_ext * read_length) / 100) ) {

                        new_header_fasta = header_fasta + "_and_" + idd2;
                        
                        if(new_header_fasta.length() > max_header_size) // //27dec 21  to avoid error in cufflink due to long ID size
                            new_header_fasta = header_fasta;
                        
                        merged_count++;

                        //merged_contigs.push_back(new_header_fasta);
                        // 28oct21 // write for chimera
                        merged_contigs.insert(make_pair(new_header_fasta, start_overlap));
                        merged_contigs.insert(make_pair(new_header_fasta, end_overlap));

                        log_new << contig_RI2 << "\t" << contig_MI2 << "\t" << "RI(FR)," << "\t" << "overlap=" << start_overlap << "-" << end_overlap << endl;
                    } else {
                        merged_check.insert(contig_RI2); // so that it won't get read again to write
                        merged_check.insert(contig_MI2);
                        if (RI_largest) { // tag RI seq, discard MI

                            string contig_RI_new = contig_RI2 + "_" + stop_RCAP3_ext;
                            mergedRI << ">" << contig_RI_new << endl << RI_seq2 << endl;

                        } else {

                            string contig_MI_new = contig_MI2 + "_" + stop_LCAP3_ext;
                            mergedRI << ">" << contig_MI_new << endl << MI_seq2 << endl;
                        }
                        merged.clear();

                    }
                } else {

                    new_header_fasta = header_fasta + "_and_" + idd2;
                    
                    if(new_header_fasta.length() > max_header_size) // //27dec 21  to avoid error in cufflink due to long ID size
                            new_header_fasta = header_fasta;
                    
                    merged_count++;

                    //merged_contigs.push_back(new_header_fasta);
                    // 28oct21 // write for chimera
                    merged_contigs.insert(make_pair(new_header_fasta, start_overlap));
                    merged_contigs.insert(make_pair(new_header_fasta, end_overlap));

                    log_new << contig_RI2 << "\t" << contig_MI2 << "\t" << "RI(FR)," << "\t" << "overlap=" << start_overlap << "-" << end_overlap << endl;
                }
            }


        } else if (!MI_dir2) { //F // take reverse complement
            // cout << "FF" << endl;
            string RI_seq_RC = utils.Rcomplement(RI_seq2); //now its R
            //MI is for means RI Is R and previously merged at left 213 so make it 312

            RI_start2 = RI_seq2.size() - (old_seq.size() - RI_start2); //9july check
            RI_end2 = RI_seq2.size() - (old_seq.size() - RI_end2);

            int ri_end = RI_seq2.size() - RI_start2;
            int ri_start = RI_seq2.size() - RI_end2;

            if (RI_seq2.size() < old_seq.size()) {
                ri_start = 0; //RI_seq.size() - RI_end;
            }

            merged = overlap_merge(MI_seq2, RI_seq_RC, MI_start2, MI_end2, ri_start, ri_end); //, MI_pos);

            if (!(merged.empty())) {


                if ((contig_RI2.find(ROASTcap3_left_tag) != std::string::npos) || (contig_MI2.find(ROASTcap3_left_tag) != std::string::npos) || (contig_RI2.find(ROASTcap3_right_tag) != std::string::npos) || (contig_MI2.find(ROASTcap3_right_tag) != std::string::npos)) {

                    int init_size;
                    bool RI_largest;

                    if (RI_seq2.size() > MI_seq2.size()) {
                        init_size = RI_seq2.size();
                        RI_largest = true;
                    } else {
                        init_size = MI_seq2.size();
                        RI_largest = false;
                    }
                    if (merged.size() > init_size && (merged.size() - init_size >= (min_CAP3_ext * read_length) / 100)) {

                        new_header_fasta = header_fasta + "_and_" + idd2;
                        
                        if(new_header_fasta.length() > max_header_size) // //27dec 21  to avoid error in cufflink due to long ID size
                            new_header_fasta = header_fasta;
                        
                        merged_count++;
                        //merged_contigs.push_back(new_header_fasta);
                        // 28oct21 // write for chimera
                        merged_contigs.insert(make_pair(new_header_fasta, start_overlap));
                        merged_contigs.insert(make_pair(new_header_fasta, end_overlap));
                        dir_update.insert(contig_RI2);

                        log_new << contig_RI2 << "\t" << contig_MI2 << "\t" << "RI(FF)," << "\t" << "overlap=" << start_overlap << "-" << end_overlap << endl;

                    } else {
                        merged_check.insert(contig_RI2); // so that it won't get read again to write
                        merged_check.insert(contig_MI2);
                         if (RI_largest) { // tag RI seq, discard MI

                            string contig_RI_new = contig_RI2 + "_" + stop_RCAP3_ext;
                            mergedRI << ">" << contig_RI_new << endl << RI_seq2 << endl;

                        } else {

                            string contig_MI_new = contig_MI2 + "_" + stop_LCAP3_ext;
                            mergedRI << ">" << contig_MI_new << endl << MI_seq2 << endl;
                        }
                        merged.clear();

                    }
                } else {

                    new_header_fasta = header_fasta + "_and_" + idd2;
                    
                     if(new_header_fasta.length() > max_header_size) // //27dec 21  to avoid error in cufflink due to long ID size
                            new_header_fasta = header_fasta;
                    
                    merged_count++;
                    //merged_contigs.push_back(new_header_fasta);
                    // 28oct21 // write for chimera
                    merged_contigs.insert(make_pair(new_header_fasta, start_overlap));
                    merged_contigs.insert(make_pair(new_header_fasta, end_overlap));
                    dir_update.insert(contig_RI2);

                    log_new << contig_RI2 << "\t" << contig_MI2 << "\t" << "RI(FF)," << "\t" << "overlap=" << start_overlap << "-" << end_overlap << endl;
                }
            }

        }


    } else if (RI_dir2) { //R

        if (MI_dir2) { //R
            //  cout << "RR" << endl;
            string RI_seq_RC = utils.Rcomplement(RI_seq2); //now its F
            // dir_update.insert(contig_RI);
            //MI is in left 3- take  RC of 12 into 21 so 213
            //  RI_start = RI_seq.size() - (old_seq.size() - RI_start);
            // RI_end = RI_seq.size() - (old_seq.size() - RI_end);

            // cout <<  RI_seq.size() << endl;
            //cout <<  MI_seq.size() << endl;
            int ri_end = RI_seq2.size() - RI_start2;
            int ri_start = RI_seq2.size() - RI_end2;

            if (RI_seq2.size() < old_seq.size()) {
                ri_start = 0; //RI_seq.size() - RI_end;
            }

            merged = overlap_merge(RI_seq_RC, MI_seq2, ri_start, ri_end, MI_start2, MI_end2); //, RI_pos);

            if (!(merged.empty())) {

                if ((contig_RI2.find(ROASTcap3_left_tag) != std::string::npos) || (contig_MI2.find(ROASTcap3_left_tag) != std::string::npos) || (contig_RI2.find(ROASTcap3_right_tag) != std::string::npos) || (contig_MI2.find(ROASTcap3_right_tag) != std::string::npos)) {

                    int init_size;
                    bool RI_largest;

                    if (RI_seq2.size() > MI_seq2.size()) {
                        init_size = RI_seq2.size();
                        RI_largest = true;
                    } else {
                        init_size = MI_seq2.size();
                        RI_largest = false;
                    }
                    if (merged.size() > init_size  && ( merged.size() - init_size >= (min_CAP3_ext * read_length) / 100)) {

                        new_header_fasta = header_fasta + "_and_" + idd2;
                        
                         if(new_header_fasta.length() > max_header_size) // //27dec 21  to avoid error in cufflink due to long ID size
                            new_header_fasta = header_fasta;
                        
                        merged_count++;
                        dir_update.insert(contig_RI2);
                         //merged_contigs.push_back(new_header_fasta);
                        // 28oct21 // write for chimera
                        merged_contigs.insert(make_pair(new_header_fasta, start_overlap));
                        merged_contigs.insert(make_pair(new_header_fasta, end_overlap));

                        log_new << contig_RI2 << "\t" << contig_MI2 << "\t" << "RI(RR)," << "\t" << "overlap=" << start_overlap << "-" << end_overlap << endl;

                    } else {
                        merged_check.insert(contig_RI2); // so that it won't get read again to write
                        merged_check.insert(contig_MI2);
                        if (RI_largest) { // tag RI seq, discard MI

                            string contig_RI_new = contig_RI2 + "_" + stop_RCAP3_ext;
                            mergedRI << ">" << contig_RI_new << endl << RI_seq2 << endl;

                        } else {

                            string contig_MI_new = contig_MI2 + "_" + stop_LCAP3_ext;
                            mergedRI << ">" << contig_MI_new << endl << MI_seq2 << endl;
                        }
                        merged.clear();

                    }
                } else {
                    new_header_fasta = header_fasta + "_and_" + idd2;
                    
                     if(new_header_fasta.length() > max_header_size) // //27dec 21  to avoid error in cufflink due to long ID size
                            new_header_fasta = header_fasta;
                    
                    merged_count++;
                    dir_update.insert(contig_RI2);
                    //merged_contigs.push_back(new_header_fasta);
                    // 28oct21 // write for chimera
                    merged_contigs.insert(make_pair(new_header_fasta, start_overlap));
                    merged_contigs.insert(make_pair(new_header_fasta, end_overlap));

                    log_new << contig_RI2 << "\t" << contig_MI2 << "\t" << "RI(RR)," << "\t" << "overlap=" << start_overlap << "-" << end_overlap << endl;
                }
            }

        } else if (!MI_dir2) { //F
          //  cout << "RF" << endl;
            // RI_start = RI_seq.size() - (old_seq.size() - RI_start);
            // RI_end = RI_seq.size() - (old_seq.size() - RI_end);
            //cout << "RI size " << RI_seq.size() << endl;
            //cout << "MI size" << MI_seq.size() << endl;
            merged = overlap_merge(MI_seq2, RI_seq2, MI_start2, MI_end2, RI_start2, RI_end2); //, MI_pos);

            if (!(merged.empty())) {
                if ((contig_RI2.find(ROASTcap3_left_tag) != std::string::npos) || (contig_MI2.find(ROASTcap3_left_tag) != std::string::npos) || (contig_RI2.find(ROASTcap3_right_tag) != std::string::npos) || (contig_MI2.find(ROASTcap3_right_tag) != std::string::npos)) {

                    int init_size;
                    bool RI_largest;

                    if (RI_seq2.size() > MI_seq2.size()) {
                        init_size = RI_seq2.size();
                        RI_largest = true;

                    } else {
                        init_size = MI_seq2.size();
                        RI_largest = false;
                    }
                    if (merged.size() > init_size  && ( merged.size() - init_size >= (min_CAP3_ext * read_length) / 100)) {

                        new_header_fasta = header_fasta + "_and_" + idd2;
                        
                         if(new_header_fasta.length() > max_header_size) // //27dec 21  to avoid error in cufflink due to long ID size
                            new_header_fasta = header_fasta;
                        
                        merged_count++;
                        //merged_contigs.push_back(new_header_fasta);
                        // 28oct21 // write for chimera
                        merged_contigs.insert(make_pair(new_header_fasta, start_overlap));
                        merged_contigs.insert(make_pair(new_header_fasta, end_overlap));

                    log_new << contig_RI2 << "\t" << contig_MI2 << "\t" << "RI(RF)," << "\t" << "overlap=" << start_overlap << "-" << end_overlap << endl;
                    } else {
                        merged_check.insert(contig_RI2); // so that it won't get read again to write
                        merged_check.insert(contig_MI2);
                         if (RI_largest) { // tag RI seq, discard MI

                            string contig_RI_new = contig_RI2 + "_" + stop_RCAP3_ext;
                            mergedRI << ">" << contig_RI_new << endl << RI_seq2 << endl;

                        } else {

                            string contig_MI_new = contig_MI2 + "_" + stop_LCAP3_ext;
                            mergedRI << ">" << contig_MI_new << endl << MI_seq2 << endl;
                        }
                        merged.clear();

                    }
                } else {
                    new_header_fasta = header_fasta + "_and_" + idd2;
                    
                     if(new_header_fasta.length() > max_header_size) // //27dec 21  to avoid error in cufflink due to long ID size
                            new_header_fasta = header_fasta;
                    
                    merged_count++;
                    //merged_contigs.push_back(new_header_fasta);
                    // 28oct21 // write for chimera
                    merged_contigs.insert(make_pair(new_header_fasta, start_overlap));
                    merged_contigs.insert(make_pair(new_header_fasta, end_overlap));

                    log_new << contig_RI2 << "\t" << contig_MI2 << "\t" << "RI(RF)," << "\t" << "overlap=" << start_overlap << "-" << end_overlap << endl;
                }
            }

        }

    }

    /*******updating old values********/
   // cout << "before updating values" << endl;
    if (!merged.empty()) {
        merged_fragments.insert(make_pair(contig_MI2, new_header_fasta)); // add new entry of contig_MI into merged_fragments

        vector<string> fragments_to_delete;
        vector<string> temp = it_newid->second;

        for (map<string, string>::iterator it_ri2 = merged_fragments.begin(); it_ri2 != merged_fragments.end(); it_ri2++) { //find all contigs present in header fasta id map and update them in merged_fragments
            vector<string >::iterator it_v = find(temp.begin(), temp.end(), it_ri2->first);

            if (it_v != temp.end()) {
                merged_fragments.insert(make_pair(it_ri2->first, header_fasta)); /// add updated header and erase old one
                fragments_to_delete.push_back(it_ri2->first); // save the id for deletion
                //merged_fragments.erase(it_ri2);

            }
        }
        //cout << "after updating merged_fragments" << endl;

        //delete the fragments marked for deletion
        for (int i = 0; i < fragments_to_delete.size(); i++) {
            //map<string, string>::iterator it_del = merged_fragments.find(fragments_to_delete[i]);
            merged_fragments.erase(fragments_to_delete[i]);

        }
        //cout << "removed old values from merged_fragments" << endl;

        //update id map with new id and add new entry in vector and erase the old one 
        temp.push_back(contig_MI2);
        id_map.insert(make_pair(new_header_fasta, temp));
        id_map.erase(it_newid);
       // cout << "updated and erased id_map" << endl;

        // update final_merged_frags both new ids and updated seq and erase the old one
        final_merged_frags.insert(make_pair(new_header_fasta, merged));
        final_merged_frags.erase(it_update_seq);
        // cout << "updated and erased final_merged_frags" << endl;

        merged_check.insert(contig_RI2);
        merged_check.insert(contig_MI2); // to check for writing in final fasta file along with other unmerged contigs
        frgaments_merged = true;
        /*******updating old values********/

        if (merged.size() <= RI_seq2.size()) // to check if size of contig didn't increase after merging of fragmented contigs - > embeded sequences
            embed_count++;

    }

    return merged_count;

}

string byRI::overlap_merge(string RI_seq, string MI_seq, int RI_start, int RI_end, int MI_start, int MI_end) {
    //cout << "overlap merge start" << endl;
    utils utils;
    int offset = (read_length*5)/100, left_boundary, right_boundary;
    
    string delimiter = "\t";
    string query_file = path_inter + "/query.fa";
    string subject_file = path_inter + "/subject.fa";
    string blast_output = path_inter + "/blast_overlap";
    string final = "";
    string ri_seq, mi_seq;

   // cout << "RI size " << RI_seq.size() << endl;
   // cout << "MI size" << MI_seq.size() << endl;

   if (RI_end > RI_seq.size()) // to avoid out of bound exception in the case of softclipping
         RI_end = RI_seq.size();
   
     if (MI_end > MI_seq.size()) // to avoid out of bound exception in the case of softclipping
         MI_end = MI_seq.size();

    if (contig_boundary) { // user def parameter to decide boundary check for overlap

        left_boundary = RI_seq.size();
        right_boundary = 0;
    } else { // consider read island boundary
        left_boundary = RI_end;
        right_boundary = MI_start;
    }

    ri_seq = RI_seq.substr(RI_start, RI_seq.size() - RI_start);

    mi_seq = MI_seq.substr(0, MI_end); //(MI_start, MI_start - MI_end);

    ofstream fasta_subject, fasta_query;

    fasta_subject.open(subject_file.c_str());
    fasta_subject << ">target" << endl << mi_seq << endl;
    fasta_subject.close();

    fasta_query.open(query_file.c_str());
    fasta_query << ">query" << endl << ri_seq << endl;
    fasta_query.close();

  //  string overlap_byBLAT = exe_path + "external_tools/blast -stepSize=5 -repMatch=2253 -minScore=5 -minIdentity=5 " + subject_file + " " + query_file + " " + blast_output + "  > /dev/null 2>&1";
   // std::system(overlap_byBLAT.c_str());
    // blastn  -evalue 5  -word_size 4 -outfmt 7 qseqid qlen sseqid slen qcovs length mismatch qstart qend sstart send' -subject -query -out 
    
    string overlap_byBLAST = "blastn -evalue 1 -word_size 4 -outfmt '7 qseqid qlen sseqid slen qcovhsp nident length mismatch qstart qend sstart send' -subject " + subject_file + " -query " + query_file + " -out " + blast_output + "  > /dev/null 2>&1";
    std::system(overlap_byBLAST.c_str());
    
    string entry;
    vector<string> temp;

    ifstream blast_left(blast_output.c_str());

    while (!(getline(blast_left, entry).eof())) {

        if (entry[0] == '#') {
            continue;
        } else { // if hit found, check score and ID of each hit
            //temp = utils.split(entry2, "\t");
            temp.clear();
            utils.str_split(entry, temp, delimiter);

            int target_start = atoi(temp[10].c_str());
            int target_end = atoi(temp[11].c_str());
            int query_size = ri_seq.size();
            int target_aln_length = target_end - target_start;
            
          //  string strand = temp[8];
            int query_start = atoi(temp[8].c_str());
            int query_end = atoi(temp[9].c_str());
 	    int q_end_contig = RI_start + query_end ;// provided query seq end and its position on contig is different

            int query_aln_length = query_end - query_start;
                       
            int score = atoi(temp[5].c_str()); // matches
            int final_score;
            
            if(query_aln_length < target_aln_length)
                final_score = float(score / query_aln_length) * 100; // score1 is number of matches / query size which is query coverage
            else
                final_score = float(score / target_aln_length) * 100;
            
	    //if (score >= 10 &&  final_score >= blast_score_TH && MI_start + target_start <= MI_start  + offset && RI_start + query_end >= RI_end - offset) { //query_size - 30) {
           // if (score >= 10 && final_score >= blast_score_TH && target_start <= 0 + offset && q_end_contig >= RI_seq.size() - offset) { //query_size - 30) { // 27/06/21 

            if (score >= min_overlap_TH && final_score >= blast_score_TH && target_start <= right_boundary + offset && q_end_contig >= left_boundary - offset) { //8th july21

                final = RI_seq.substr(0, q_end_contig) + MI_seq.substr(target_end, MI_seq.size() - target_end);
                
                if (!final.empty()) // update merged contigs map with start and end positions of overlap
                {
                    //28oct21 overlap boundary extraction to avoid split of same position in multi-gene chimera
                    start_overlap = RI_start; // start of overlap
                    end_overlap = q_end_contig; // end of overlap

                }
            }
            break;
        } 

    }

    blast_left.close();
    remove(query_file.c_str());
    remove(blast_output.c_str());
    remove(subject_file.c_str());
    //cout << "overlap merge end" << endl;


    return final;
}

/*string byRI::overlap_merge(string str1, string str2) {//, int str1_pos) {
    //cout << "overlap merge start" << endl;
    utils utils;
    int discard_seq_TH = 20;
    int i = str1.size() - read_length, n, m, pos = 1, pos1;
    int j, max_mismatch = 10; //25 commented on 10 feb
    int k = 1;
    int counter = 0, mismatch = 0;
    string overlap, maxOverlap, final, MI, RI, truncated_RI, truncated_MI;
    int overlapSize = 0;
    // cout << str1.size() << " " << str2.size() << endl;
    
    for (; i < str1.size(); i++) {
        for (j = 0; j <= read_length + 25; j++) { //str2.size() / 2
            for (k = 1; k <= read_length ; k++) { //str2.size()
                if (str1.substr(i, k) == str2.substr(j, k)) {
                    overlap = str1.substr(i, k);

                } else {
                    if (overlap.size() > maxOverlap.size()) {
                        maxOverlap = overlap;
                        pos = j + k;
                        pos1 = i;
                        RI = str1.substr(0, pos1);
                        MI = str2.substr(pos - 1, str2.size() - 1);
                        truncated_RI = str1.substr(pos1 + maxOverlap.size(), str1.size() - 1);
                        truncated_MI = str2.substr(0, (pos - maxOverlap.size()));
                        //std::transform(maxOverlap.begin(), maxOverlap.end(), maxOverlap.begin(), ::tolower);
                        final = RI + maxOverlap + MI;
                        overlap.clear();
                    }
                }
            }
        }
    }
    char base_str1, base_str2;
    int same = 0; // to check consecutive 3 mismatches
    if (final != "") {
        if (truncated_RI.size() > discard_seq_TH) {
            //cout << "complete truncated_RI " << endl;
            i = pos1 - 1; // + maxOverlap.size();
            j = pos - 1;
            
            while (i < str1.size()) {
                
                while (j < str2.size()) {
                    base_str1 = str1[i];
                    base_str2 = str2[j];
                    if (base_str1 == base_str2) {
                        overlap = overlap + base_str1;
                        same = 0;
                        i++;
                        j++;
                        goto out_of_j;
                    } else {

                        char code = utils.base_code(base_str1, base_str2);
                        mismatch++;
                        
                        if (mismatch <= max_mismatch && same <= 5) {
                            overlap = overlap + code;
                            i++;
                            j++;
                            pos = j;
                            same++;
                            goto out_of_j;
                        } else {
                            maxOverlap = maxOverlap + overlap;
                            pos = j; // size of overlap
                            pos1 = i;
                            RI = str1.substr(0, pos1);
                            MI = str2.substr(pos - 1, str2.size() - 1);
                            truncated_RI = str1.substr(pos1, str1.size() - 1);
                            truncated_MI = str2.substr(0, (pos - overlap.size()));
                            // std::transform(maxOverlap.begin(), maxOverlap.end(), maxOverlap.begin(), ::tolower);
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
        } else if (truncated_MI.size() > discard_seq_TH) {
            // cout << "complete truncated_MI " << endl;
            i = pos1 - maxOverlap.size();
            j = pos - maxOverlap.size();
            while (i >= 0) {
                while (j >= 0) {
                    if (j == 0) goto out_of_i;
                    base_str1 = str1[i];
                    base_str2 = str2[j];
                    if (base_str1 == base_str2) {
                        overlap = base_str1 + overlap;
                        same = 0;
                        i--;
                        j--;
                        goto out_of_jj;
                    } else {

                        char code = utils.base_code(base_str1, base_str2);
                        mismatch++;
                        if (mismatch <= max_mismatch && same <= 5) {
                            overlap = code + overlap;
                            i--;
                            j--;
                            pos = j;
                            same++;
                            goto out_of_jj;
                        } else {
                            maxOverlap = overlap + maxOverlap;
                            pos = (j - 1) + (maxOverlap.size() - 1); //str1.size() - maxOverlap.size();
                            pos1 = i - 1;
                            RI = str1.substr(0, pos1);
                            MI = str2.substr(pos, str2.size() - 1);
                            truncated_RI = str1.substr(pos1, str1.size() - 1);
                            truncated_MI = str2.substr(0, (pos - overlap.size()));
                            // std::transform(maxOverlap.begin(), maxOverlap.end(), maxOverlap.begin(), ::tolower);
                            final = RI + maxOverlap + MI;
                            overlap.clear();
                            goto out_of_ii;
                        }
                    }
                }

out_of_jj:
                ;
            }
out_of_ii:
            ;


        } else {// ((truncated_RI.size() < discard_seq_TH) && (truncated_MI.size() < discard_seq_TH)) {
            //cout << "truncated RI:" << truncated_RI << endl;
            //cout << "truncated MI:" << truncated_MI << endl;
            return final;
        }
        
        if ((truncated_RI.size() < discard_seq_TH) && (truncated_MI.size() < discard_seq_TH)) {
            // cout << "truncated RI:" << truncated_RI << endl;
            // cout << "truncated MI:" << truncated_MI << endl;
            return final;
        } else {
            final = "";
            return final;
        }
    } else if (final == "" && overlap != "") {
        maxOverlap = overlap;
        RI = str1.substr(0, i);
        MI = str2.substr(j - 1, str2.size() - 1);
        //std::transform(maxOverlap.begin(), maxOverlap.end(), maxOverlap.begin(), ::tolower);
        final = RI + maxOverlap + MI;
        overlap.clear();
        return final;

    } else {
        // cout << "no overlap found" << endl;
        //cout << " overlap merge end" << endl;

        return final;
    }
}*/

void byRI::extend_unmappedBam(string bam_file, string in_file, string &out_file, string &log) {

    //typedef vector<vector<string> > all_vectors;
    // all_vectors all_vec;
    int contig_count = 0, count = 1, ext_count = 0;

    if (unmapped_contig.size() > 0) { 
        int THREADS = atoi(threads.c_str());
        int num_file_generated = 1;

        vector<string> temp_vec;
        std::stringstream ct;
        ct << count;

        string file_vec = path_inter + "/vect_file_" + ct.str() + ".fasta"; // write for each 
        ct.str(string());

        ofstream vec_file;
        vec_file.open(file_vec.c_str());

      //  cout << unmapped_contig.size() << endl;
        int sub_vec_size;

        if (unmapped_contig.size() > THREADS && THREADS > 1)
            sub_vec_size = unmapped_contig.size() / THREADS;
        else
            sub_vec_size = unmapped_contig.size();

     //   cout << sub_vec_size << endl;
        //int remaining_contigs = contig_count % allowed_threads; // to deal with odd number of contigs

        //if(remaining_contigs > 0)
        //    sub_vec_size = sub_vec_size + 1;

        //write contigs names in n number of files based on number of threads
        vector<string>::iterator i = unmapped_contig.begin();

        while (i != unmapped_contig.end()) {

            if (count <= THREADS - 1) {

                if (contig_count < sub_vec_size) { // write n contig ids in current file

                    vec_file << *i << endl;
                    //temp_vec.push_back(*i);
                    contig_count++;

                } else { // close current file and open new file for n contigs again
                    vec_file.close();
                    count++;
                    ct << count;

                    string file_vec = path_inter + "/vect_file_" + ct.str() + ".fasta"; // write for each 
                    vec_file.open(file_vec.c_str());
                    num_file_generated++;

                    //  all_vec.push_back(temp_vec);
                    //temp_vec.clear();
                    contig_count = 0;
                    ct.str(string());

                }
            } else { //add all remaining data in last file
                vec_file << *i << endl;
                //temp_vec.push_back(*i);
            }
            i++;
        }
        unmapped_contig.clear();
        vec_file.close();
        ct.str(string());

        string ROAST_extendContigs;
        stringstream P1, P2, P3, P4;

        //cout << "before exe call" << endl;
        P1 << min_overlap_TH;
        P2 << min_unmapped_reads_CAP3;
	P3 << min_CAP3_ext;
        P4 << read_length;
       // cout << ct.str() << endl;
        //call TRAT_extendCOntigs.exe for all input IDs files generated previously and send in the backgroud to run
        for (int i = 1; i <= num_file_generated; i++) {
            ct << i;
            //  generate cap3 assemblies for contigs having unmapped reads of terminally mapped mates
            // ROAST_extendContigs tt;
            // tt.cap3_extension(bam_file, in_file, log, ct.str(), path_inter, exe_path, min_overlap_TH, min_unmapped_reads_CAP3, min_CAP3_ext, read_length, tool_name_flag);

            ROAST_extendContigs = exe_path + "/ROAST_extendContigs " + bam_file + " " + in_file + " " + log + " " + ct.str() + " " + path_inter + " " + exe_path + " " + P1.str() + " " + P2.str() + " " + P3.str() +  " " + P4.str() +" " +tool_name_flag + " &";
            std::system(ROAST_extendContigs.c_str());

           ct.str(string());
        }
        P1.str(string());
        P2.str(string());
	P3.str(string());
        P4.str(string());
        ct.str(string());
        //  cout << "after exe call before deletion" << endl;

        // at the end of RAT_extendCOntigs write confirmation text file to ensure job in that thread has been executed, later check status of those confirmation files to move to next step
        int running_threads_count = num_file_generated;
        //cout << ct.str() << endl;
        while (running_threads_count > 0) {

            for (int i = 1; i <= num_file_generated; i++) {
                ct << i;
                string thread_done = path_inter + "/done" + ct.str() + ".txt";
                 ct.str(string());

                if (boost::filesystem::exists(thread_done)) // does filePath actually exist?
                {
                    running_threads_count--;
                } else {
                    //i = 1; // start from 1st file again to check status
                    sleep(10); // 
                    running_threads_count = num_file_generated;
                    break;
                }

            }
        }
        ct.str(string());
        //  cout << " all threads done" << endl;

        //merge all new_assemblies


        string merge_assemblies = "cat ";

        for (int i = 1; i <= num_file_generated; i++) {
            ct << i;
            string individual_extended_files = path_inter + "/ROAST_extendContigs_" + ct.str() + ".fasta";
            merge_assemblies = merge_assemblies + individual_extended_files + " ";
            ct.str(string());
        }
        merge_assemblies = merge_assemblies + " > " + out_file;

        std::system(merge_assemblies.c_str());
        ct.str(string());
        // cout << "sub assemblies merged " << endl;

        //remove all intermediate file of this process
        for (int i = 1; i <= num_file_generated; i++) {

            ct << i;
            string thread_done = path_inter + "/done" + ct.str() + ".txt";
            string individual_extended_files = path_inter + "/ROAST_extendContigs_" + ct.str() + ".fasta";
            remove(thread_done.c_str());
            remove(individual_extended_files.c_str());
            ct.str(string());

        }

        //  cout << "fasta and done files removed" << endl;
        ct.str(string());
        string entry_ID;
        for (int i = 1; i <= num_file_generated; i++) {
            //cap3_extension also adds contigs_name in cap3_left_list, cap3_right_list for unmerged assemblies and merged_contigs if has been merged with main contig
            ct << i;
            string extended_byCAP3_file = path_inter + "/extended_byCAP3_" + ct.str() + ".txt";
            ifstream ext;
            ext.open(extended_byCAP3_file.c_str());
            bool get_count = true;

            while (!(getline(ext, entry_ID).eof())) {

                if (get_count) { // first line of extended_byCAP3 file contains number of total extensions 
                    ext_count = ext_count + atoi(entry_ID.c_str());
                    get_count = false;
                }
                extended_byCAP3.insert(entry_ID);

            }
            ext.close();
            remove(extended_byCAP3_file.c_str());

            string merged_contigs_file = path_inter + "/merged_contigs_" + ct.str() + ".txt";
            ifstream merged_file;
            merged_file.open(merged_contigs_file.c_str());

            /*while (!(getline(merged_file, entry_ID).eof())) {

                merged_contigs.push_back(entry_ID);

            }*///28oct21
            merged_file.close();
            remove(merged_contigs_file .c_str());

            string cap3_left_file = path_inter + "/cap3_left_" + ct.str() + ".txt";
            ifstream left_file;
            left_file.open(cap3_left_file.c_str());

            while (!(getline(left_file, entry_ID).eof())) {

                cap3_left_list.push_back(entry_ID);

            }
            left_file.close();
            remove(cap3_left_file.c_str());

            string cap3_right_file = path_inter + "/cap3_right_" + ct.str() + ".txt";
            ifstream right_file;
            right_file.open(cap3_right_file.c_str());


            while (!(getline(right_file, entry_ID).eof())) {

                cap3_right_list.push_back(entry_ID);

            }

            right_file.close();
            remove(cap3_right_file.c_str());

            //set <string> extended_byCAP3;    list<string> merged_contigs, cap3_left_list, cap3_right_list;
            string file_vec = path_inter + "/vect_file_" + ct.str() + ".fasta"; // write for each 
            remove(file_vec.c_str());
            ct.str(string());
        }
        ct.str(string());
    } // if unmapped_contigs > 0

    cout << "No. of extended contigs by Cap3 at iteration " << ": " << ext_count << endl;
}
































