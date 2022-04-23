/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ROAST_extendContigs.cpp
 * Author: madiha
 * 
 * Created on December 23, 2020, 4:03 PM
 */

#include "ROAST_extendContigs.h"
//#include "global.h"
#include "utils.h"
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
//#include <boost/regex.hpp>

ROAST_extendContigs::ROAST_extendContigs() {
}

ROAST_extendContigs::ROAST_extendContigs(const ROAST_extendContigs& orig) {
}

ROAST_extendContigs::~ROAST_extendContigs() {
}

string ROAST_LE_tag = "trLE";

string ROAST_RE_tag = "trRE";

string ROAST_BE_tag = "trBE";

string ROASTcap3_left_tag = "CAP3LT";

string ROASTcap3_right_tag = "CAP3RT";

string stop_LCAP3_ext = "CAP3LX";

string stop_RCAP3_ext = "CAP3RX";

void ROAST_extendContigs::cap3_extension(string bam_file, string in_file,  string &log, string num, string path_inter, string exe_path, int min_overlap_TH, int min_unmappead_reads_CAP3, int min_CAP3_ext, int read_length, string tool_name_flag) { 

    utils utils;
    //std::size_t found = bam_file.find_last_of("/\\");
    bool merged = false, left = false, right = false;
    string entry_ID;
    
    string file_vec = path_inter + "/vect_file_" + num + ".fasta"; // write for each 
    ifstream vec_file;
    vec_file.open(file_vec.c_str());
    
    string unmapped_reads = path_inter+ "/unmapped_reads." + num + ".bam";
    string unmapped_mates = path_inter + "/unmapped_mates." + num + ".bam";
    string debugFile = path_inter + "/debug_cap3_MiassemblyIDs" + num + ".txt";
    
    string fa_file = path_inter + "/unmapped1.fa" + num;
    string out_file = path_inter + "/ROAST_extendContigs_" + num + ".fasta"; // write for each 
    string debug_file = path_inter + "/debug_ROAST_extendContigs_" + num + ".txt";
    string assembly = fa_file + ".cap.contigs";
    
   // ofstream debug;
   // debug.open(debug_file.c_str());

   //string ext_assem = bam_file.substr(0, found) + "/extended.cap3_1.check_newFormula.fasta";


   // file to get data for global varaibles set <string> extended_byCAP3; list<string> merged_contigs, cap3_left_list, cap3_right_list;

    set <string> extended_byCAP3;
    list<string> merged_contigs, cap3_left_list, cap3_right_list;
     
    ofstream extended_assembly;
    extended_assembly.open(out_file.c_str());
  
    ofstream log_cap3;
    log_cap3.open(log.c_str(), ios::app);

    int prv_ref_id = 0, cur_ref_id, unmapped_seq = 0, temp_contigs = 0, ext_count = 0;
    string contig_name, entry, left_new_contig, right_new_contig, original_contig, new_assembly, flag, cap3_header;

   // string left_cap3, right_cap3;
    int len_read_contig, contig_length;

    /*extract unmapped reads from each contig one by one
       make assembly from cap3 
       find overlap with existing contig to extend it*/

        while (!(getline(vec_file, contig_name).eof())) {             

        if ((contig_name.find(stop_LCAP3_ext) == std::string::npos)) { // contig doesnt have tag to stop making left CAP3 (to avoid infinite loop)
           
            left_new_contig = left_cap3_assembly(bam_file, debugFile, unmapped_reads, unmapped_mates, fa_file, assembly, contig_name, exe_path, min_overlap_TH, min_unmappead_reads_CAP3);
        }


        original_contig = utils.extract_fasta(contig_name, in_file);

        int len = original_contig.length();
        stringstream contig_len;
        contig_len << len;

        if ((contig_name.find(stop_RCAP3_ext) == std::string::npos)) { // contig doesnt have tag to stop making right CAP3 (to avoid infinite loop)
           
            right_new_contig = right_cap3_assembly(bam_file, debugFile, unmapped_reads, unmapped_mates, fa_file, assembly, contig_name, exe_path, contig_len.str(), min_overlap_TH, min_unmappead_reads_CAP3);
        }

        //left_new_contig = left_cap3; //it_left->first;
        //right_new_contig = right_cap3; //it_right->first;

        size_t tag = contig_name.find(tool_name_flag);

        if (tag != string::npos) // tag found
            cap3_header = contig_name.substr(0, tag - 1); // remove tag
        else
            cap3_header = contig_name;

        if (!(left_new_contig.empty()) && !(right_new_contig.empty())) {
          //  cout << contig_name << endl;
            
            // Ns = it_left->second;
            flag = "left";
            new_assembly = overlap_mergeCAP3(left_new_contig, original_contig, flag, min_overlap_TH, path_inter, exe_path, num);

            if (new_assembly.empty()) { // no overlap found between new contig and existing contig

                merged = false;
                extended_assembly << ">" << cap3_header << "_" << ROASTcap3_left_tag << endl << left_new_contig << endl;
                log_cap3 << contig_name << "\t" << "cap3 left_side_notMerged" << endl;
                cap3_left_list.push_back(cap3_header);
                new_assembly = original_contig;
                temp_contigs++;
                extended_byCAP3.insert(contig_name);

            } else {
                // if new_assembly size is more than 50% of read length to the initial contig consider overlap, otherwise add stopCAP3extension tag to the contig and drop extension
                if (new_assembly.size() > original_contig.size() && (new_assembly.size() - original_contig.size() >= (min_CAP3_ext * read_length) / 100) ) {
                    merged = true;
                    extended_byCAP3.insert(contig_name);
                } else {
                      extended_byCAP3.insert(contig_name);
                    contig_name = contig_name + "_" + stop_LCAP3_ext;
                    merged = false;
                    new_assembly = original_contig;
                }

            }
            flag = "right";

            string right_extension = overlap_mergeCAP3(right_new_contig, new_assembly, flag, min_overlap_TH, path_inter, exe_path, num);  

            if (right_extension.empty()) {
                
                cap3_right_list.push_back(cap3_header);
                extended_assembly << ">" << cap3_header << "_" << ROASTcap3_right_tag << endl << right_new_contig << endl;
                log_cap3 << contig_name << " cap3 right_sides_notMerged" << endl;
                temp_contigs++;

                if (!merged)
                    extended_assembly << ">" << contig_name << endl << new_assembly << endl;
                else // merged true for left but not for right  /*updated 3 nov 2020*/
                {
                    left = true; //extended_assembly << ">" << contig_name << endl << new_assembly << endl;
                    merged_contigs.push_back(contig_name);
                }
            } else { // if right ext is not empty and meet length criteria merge them, otherwise update previously updated (or not) contig ID for no ext and add in assembly

                // if new_assembly size is more than 50% of read length to the initial contig consider overlap, otherwise add stopCAP3extension tag to the contig and drop extension
                if (right_extension.size() > new_assembly.size() && (right_extension.size() - new_assembly.size() >= (min_CAP3_ext * read_length) / 100)) {

                    new_assembly = right_extension;
                    // extended_assembly << ">" << contig_name << endl << right_extension << endl;
                    merged_contigs.push_back(contig_name);
                    log_cap3 << contig_name << "\t" << "cap3 both_sides_extended" << endl;
                    right = true, left = true;

                } else {
                    extended_byCAP3.insert(contig_name);
                    contig_name = contig_name + "_" + stop_RCAP3_ext;
                    extended_assembly << ">" << contig_name << endl << new_assembly << endl;

                }


            }
            extended_byCAP3.insert(contig_name);

        }// only left side is assembled 
        else if (!(left_new_contig.empty())) {
            //  cout << contig_name << endl;
            flag = "left";
            //  Ns = it_left->second;
            new_assembly = overlap_mergeCAP3(left_new_contig, original_contig, flag, min_overlap_TH, path_inter, exe_path, num);

            if (new_assembly.empty()) { // no overlap found between new contig and existing contig

                new_assembly = original_contig;
                cap3_left_list.push_back(cap3_header);
                extended_assembly << ">" << contig_name << endl << original_contig << endl;
                extended_assembly << ">" << cap3_header << "_" << ROASTcap3_left_tag << endl << left_new_contig << endl;
                temp_contigs++;

                log_cap3 << contig_name << "\t" << "cap3 left_side_notMerged" << endl;
                // new_assembly = left_new_contig + string(Ns, 'N') + original_contig;
            } else {

                if (new_assembly.size() > original_contig.size() && (new_assembly.size() - original_contig.size() >= (min_CAP3_ext * read_length) / 100)) {
                    // extended_assembly << ">" << contig_name << endl << new_assembly << endl;

                    log_cap3 << contig_name << "\t" << "cap3 left_side_extended" << endl;
                    left = true;

                } else {
                      extended_byCAP3.insert(contig_name);
                    contig_name = contig_name + "_" + stop_LCAP3_ext;
                    extended_assembly << ">" << contig_name << endl << new_assembly << endl;

                }

            }
            extended_byCAP3.insert(contig_name);

        }// only right side is assembled 
        else if (!(right_new_contig.empty())) {
          //  cout << contig_name << endl;
            flag = "right";
            //  Ns = it_right->second;
            new_assembly = overlap_mergeCAP3(right_new_contig, original_contig, flag, min_overlap_TH, path_inter, exe_path, num);

            if (new_assembly.empty()) {
                new_assembly = original_contig;
                cap3_right_list.push_back(cap3_header);

                extended_assembly << ">" << contig_name << endl << original_contig << endl;
                extended_assembly << ">" << cap3_header << "_" << ROASTcap3_right_tag << endl << right_new_contig << endl;
                temp_contigs++;
                log_cap3 << contig_name << "\t" << "cap3 right_side_notMerged" << endl;
                // new_assembly = original_contig + string(Ns, 'N') + right_new_contig;
            } else {
                //extended_assembly << ">" << contig_name << endl << new_assembly << endl;
                if (new_assembly.size() > original_contig.size()  && (new_assembly.size() - original_contig.size() >= (min_CAP3_ext * read_length) / 100)) {

                    merged_contigs.push_back(contig_name);
                    log_cap3 << contig_name << "\t" << "cap3 right_side_extended" << endl;
                    right = true;

                } else {
                      extended_byCAP3.insert(contig_name);
                    contig_name = contig_name + "_" + stop_RCAP3_ext;
                    extended_assembly << ">" << contig_name << endl << new_assembly << endl;
                }
            }
            extended_byCAP3.insert(contig_name);
        }

        if (left || right) { //if extended

            size_t update_ind = contig_name.find(tool_name_flag);
            // size_t update_ind2 = contig_name.find_last_of(tool_name_flag);

            if (update_ind != string::npos) //[] found
            {

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
                        ext_count++; //not already extended -> to get unique count of each extended contig

                    }
                } else if (left) {

                    if (find_BE != string::npos) {
                        //ignore
                    } else if (find_RE != string::npos) {

                        contig_name = contig_name.replace(find_RE, ROAST_BE_tag.size(), ROAST_BE_tag);

                    } else if (find_LE != string::npos) {
                        //ignore
                    } else {
                        contig_name = contig_name.substr(0, contig_name.size() - tool_name_flag.length() - 1) + "-" + ROAST_LE_tag  + "_" + tool_name_flag;
                        ext_count++; //not already extended -> to get unique count of each extended contig
                    }
                } else{// if (right) {

                    if (find_BE != string::npos) {
                        //ignore
                    } else if (find_RE != string::npos) {
                        //ignore
                    } else if (find_LE != string::npos) {

                        contig_name = contig_name.replace(find_LE, ROAST_BE_tag.size(), ROAST_BE_tag);

                    } else {

                        contig_name = contig_name.substr(0, contig_name.size() - tool_name_flag.length() - 1) + "-" + ROAST_RE_tag + "_" + tool_name_flag;
                        ext_count++; //not already extended -> to get unique count of each extended contig

                    }
                }
            } else { // NO improvement yet                 
                    if (left && right)
                        contig_name = contig_name + "_" + tool_name_flag + "_" + ROAST_BE_tag+"_" + tool_name_flag;
                    else if (left)
                        contig_name = contig_name + "_" + tool_name_flag + "_" + ROAST_LE_tag +"_" + tool_name_flag;
                    else// (right == true)
                        contig_name = contig_name + "_" + tool_name_flag + "_" + ROAST_RE_tag + "_" + tool_name_flag;
                    ext_count++;
            }
            extended_assembly << ">" << contig_name << endl << new_assembly << endl;
        }
       // debug << contig_name << endl;
        
        new_assembly = "";
        left = false;
        right = false;
        remove(unmapped_reads.c_str());
        remove(fa_file.c_str());
        contig_len.str(string());
    }
    extended_assembly.close();
    log_cap3.close();
    vec_file.close();

   // extended_incomplete_contigs = extended_incomplete_contigs + ext_count;

   // cout << "No. of extended contigs by Cap3 at iteration " << ": " << ext_count << endl;
  //  cout << "Temporary contigs:" << temp_contigs << endl;

    string extended_byCAP3_file = path_inter + "/extended_byCAP3_" + num + ".txt";
    ofstream ext;
    ext.open(extended_byCAP3_file.c_str());
    
    ext << ext_count << endl;
    
    for (set<string>::iterator i = extended_byCAP3.begin(); i != extended_byCAP3.end(); ++i) {
        ext << *i << endl;
    }
    
    ext.close();
  //  debug << " extended_byCAP3 writing done " << endl;
    
    string merged_contigs_file = path_inter + "/merged_contigs_" + num + ".txt";
    ofstream merged_file;
    merged_file.open(merged_contigs_file.c_str());

    for (list<string>::iterator i = merged_contigs.begin(); i != merged_contigs.end(); i++) {
        merged_file << *i << endl;
    }
    merged_file.close();
   /// debug << " merged_file writing done " << endl;
            
    string cap3_left_file = path_inter + "/cap3_left_" + num + ".txt";
    ofstream left_file;
    left_file.open(cap3_left_file.c_str());
    
    for (list<string>::iterator i = cap3_left_list.begin(); i != cap3_left_list.end(); ++i) {
        left_file << *i << endl;
    }
    
    left_file.close();
  //  debug << " cap3_left writing done " << endl;
    
    string cap3_right_file = path_inter + "/cap3_right_" + num + ".txt";
    ofstream right_file;
    right_file.open(cap3_right_file.c_str());


    for (list<string>::iterator i = cap3_right_list.begin(); i != cap3_right_list.end(); ++i) {

        right_file << *i << endl;

    }

    extended_byCAP3.clear();
    merged_contigs.clear();
    cap3_left_list.clear();
    cap3_right_list.clear();
    
    right_file.close();
   // debug << " cap3_right writing done " << endl;
   // return extended_byCAP3;
   // debug << " job done" << endl;
    
    string thread_done = path_inter + "/done" + num + ".txt";
    ofstream done;
    done.open(thread_done.c_str());
    done << "Job Finished";
    done.close();
}



string ROAST_extendContigs::left_cap3_assembly(string bam_file, string debugFile, string unmapped_reads, string unmapped_mates, string fa_file, string assembly, string contig_name, string exe_path, int min_overlap_TH, int min_unmappead_reads_CAP3) { // we don't need Ns here so return in string
    BamReader reader;
    BamTools::BamAlignment al;

    int unmapped_seq = 0,  cap3_assembly_count = 0;
    set <string> mate_ids;
   // string temp_cap3 = path_inter + "/old_cap3_assembly.fa";

    vector<CigarOp> cigar;


    ofstream fa;// debugfile;
   //// debugfile.open(debugFile.c_str());
    
    string entry, new_contig = "";
    fa.open(fa_file.c_str());

   // string extract_read_unmappedMates = "samtools view -b -f 8 -F 4 " + bam_file + " ''" + contig_name + ":1-100'' " + " -o " + unmapped_mates + "  > /dev/null 2>&1";
    
    string extract_read_unmappedMates = "samtools view -b -f 8 -F 260 " + bam_file + " ''" + contig_name + ":1-100'' " + " -o " + unmapped_mates + "  > /dev/null 2>&1";
    string extract_unmapped_reads = "samtools view -u -f 4 -F264 " + bam_file + " ''" + contig_name  + ":1-100'' " + " -o " + unmapped_reads + "  > /dev/null 2>&1";

    system(extract_unmapped_reads.c_str()); //read s which have unmapped mates
    system(extract_read_unmappedMates.c_str()); // unmapped reads

    if (!reader.Open(unmapped_mates)) {
        cerr << "Could not open BAM file." << endl;
        exit(0);
    }

    while (reader.GetNextAlignment(al)) {
        
        if(al.IsReverseStrand()){ // read should be on reverse dir for left unmapped reads
            mate_ids.insert(al.Name);
            //count++;
        }
    } 
    
    if (!reader.Open(unmapped_reads)) {
        cerr << "Could not open BAM file." << endl;
        exit(0);
    }
    
    set<string>::iterator itr = mate_ids.begin();
    while (reader.GetNextAlignment(al)) {
         
        itr  = mate_ids.find(al.Name);
            if (itr != mate_ids.end())  {
                fa << ">" << al.Name << endl << al.QueryBases << endl;
                unmapped_seq++;
            }
    }

mate_ids.clear();
    if (unmapped_seq >= min_unmappead_reads_CAP3) { //make assembly for more than 2 unmapped reads
        string cap3_assembly = exe_path + "external_tools/cap3 " + fa_file + "  > /dev/null 2>&1";
        system(cap3_assembly.c_str());
        
        string temp = fa_file + ".cap.ace";
        remove(temp.c_str());
        temp = fa_file + ".cap.contigs.links";
        remove(temp.c_str());
        temp = fa_file + ".cap.contigs.qual";
        remove(temp.c_str());
        temp = fa_file + ".cap.info";
        remove(temp.c_str());
        temp = fa_file + ".cap.singlets";
        remove(temp.c_str());
        ifstream assem;
        assem.open(assembly.c_str());

        if (assem.peek() != std::ifstream::traits_type::eof()) //not empty file
        {

            while (!(getline(assem, entry).eof())) { //  merge with orignal contig

                if (entry[0] == '>' && cap3_assembly_count < 1) { //ignore header 
                    //continue;
                    cap3_assembly_count++;
                } else if (entry[0] == '>' && cap3_assembly_count >= 1) {
                    new_contig = ""; //debugfile << contig_name << endl;
                    goto omit_assembly;
                } else
                    new_contig = new_contig + entry;
            }

        }
                
    }
omit_assembly:
        ;

        
        /*re_assembly:;
        if (assem.peek() != std::ifstream::traits_type::eof()) //not empty file
        {

            while (!(getline(assem, entry).eof())) { //  merge with orignal contig
                if (entry[0] == '>' && cap3_assembly_count < 1) //ignore header
                {
                    //continue;
                    cap3_assembly_count++;
                }
                else if(entry[0] == '>' && cap3_assembly_count >= 1){ // re run cap3 if more than 1 contigs generated in prev run
                    string mv_command = " cp " + assembly + " " + temp_cap3;
                    system(mv_command.c_str());
                    string cap3_assembly = exe_path + "external_tools/cap3 " + temp_cap3 + "  > /dev/null 2>&1";
                    system(cap3_assembly.c_str());
                    
                    ifstream assem;
                    assem.open(assembly.c_str());
                    cap3_assembly_count++;
                    goto re_assembly;
                }
                
                else {
                    new_contig = new_contig + entry;
                }
            }

        }*/

    reader.Close();
   // debugfile.close();
    remove(unmapped_reads.c_str());
    remove(unmapped_mates.c_str());
    remove(assembly.c_str());
    return new_contig;
}



string ROAST_extendContigs::right_cap3_assembly(string bam_file, string debugFile, string unmapped_reads, string unmapped_mates, string fa_file, string assembly, string contig_name, string exe_path , string contig_length, int min_overlap_TH, int min_unmappead_reads_CAP3) {
    utils utils;
    BamReader reader;
    BamTools::BamAlignment al;

    int unmapped_seq = 0, cap3_assembly_count = 0;
    set <string> mate_ids;

    vector<int> clipsize, read_pos, gen_pos;
    vector<CigarOp> cigar;

    ofstream fa;//  debugfile;
   // debugfile.open(debugFile.c_str());
    
    string entry, new_contig = "";
    fa.open(fa_file.c_str());

   int start_pos = atoi(contig_length.c_str()) - 200;
    stringstream ss1;
    ss1 << start_pos;
    
    
    string extract_read_unmappedMates = "samtools view -b -f 8 -F 260 " + bam_file + " ''" + contig_name + ":" + ss1.str() + "-" + contig_length + "'' -o " + unmapped_mates + "  > /dev/null 2>&1";
    string extract_unmapped_reads = "samtools view -u -f 4 -F264 " + bam_file + " ''" + contig_name  + ":" + ss1.str() + "-" + contig_length + "'' -o " + unmapped_reads + "  > /dev/null 2>&1";

    system(extract_unmapped_reads.c_str());
    system(extract_read_unmappedMates.c_str());

    if (!reader.Open(unmapped_mates)) {
        cerr << "Could not open BAM file." << endl;
        exit(0);
    }

    while (reader.GetNextAlignment(al)) {
        if(!(al.IsReverseStrand())){ // read should be on for dir for right unmapped reads
            mate_ids.insert(al.Name);
           // count++;
        }
        }
  
    if (!reader.Open(unmapped_reads)) {
        cerr << "Could not open BAM file." << endl;
        exit(0);
    }

    set<string>::iterator itr = mate_ids.begin();
    while (reader.GetNextAlignment(al)) {

        itr = mate_ids.find(al.Name);
        if (itr != mate_ids.end()) {
            fa << ">" << al.Name << endl << al.QueryBases << endl;
            unmapped_seq++;
        }
    }
    mate_ids.clear();
    fa.close();
    if (unmapped_seq >= min_unmappead_reads_CAP3) { //make assembly for more than 2 unmapped reads
        string cap3_assembly = exe_path + "external_tools/cap3 " + fa_file + "  > /dev/null 2>&1";
        system(cap3_assembly.c_str());
        string temp = fa_file + ".cap.ace";
        remove(temp.c_str());
        temp = fa_file + ".cap.contigs.links";
        remove(temp.c_str());
        temp = fa_file + ".cap.contigs.qual";
        remove(temp.c_str());
        temp = fa_file + ".cap.info";
        remove(temp.c_str());
        temp = fa_file + ".cap.singlets";
        remove(temp.c_str());

        ifstream assem;
        assem.open(assembly.c_str());
        if (assem.peek() != std::ifstream::traits_type::eof()) //empty file
        {

            while (!(getline(assem, entry).eof())) { //  merge with orignal contig
                if (entry[0] == '>' && cap3_assembly_count < 1) { //ignore header  //ignore header
                    //continue;
                    cap3_assembly_count++;
                } else if (entry[0] == '>' && cap3_assembly_count >= 1) { // ignore if cap3 generates more than one assembly (can be due to library/sequencing error or paralogs)
                    new_contig = ""; //debugfile << contig_name << endl;
                    goto omit_assembly;
                } else
                    new_contig = new_contig + entry;
            }

            new_contig = utils.Rcomplement(new_contig); // take reverse complement of mates or new contig
        }
    }
    omit_assembly:
        ;

    reader.Close();
    //debugfile.close();
    remove(unmapped_reads.c_str());
    remove(unmapped_mates.c_str());
    remove(assembly.c_str());

    return new_contig;
}
string ROAST_extendContigs::overlap_mergeCAP3(string query, string target, string flag, int min_overlap_TH, string path_inter, string exe_path, string num) {
    utils utils;
    string delimiter = "\t";
    string target_file = path_inter + "/target" + num + ".fasta";
    string query_file = path_inter + "/query" + num + ".fasta";
    string blat_output = path_inter + "/blat_overlap_" + num;
    string final = "";

    ofstream fasta_subject, CAP3_query;
    fasta_subject.open(target_file.c_str());
    CAP3_query.open(query_file.c_str());
    
    fasta_subject << ">subject" << endl << target << endl;
    fasta_subject.close();

    CAP3_query<< ">CAP3_query" << endl << query << endl;
    CAP3_query.close();
    
    // string overlap_byBLAT = exe_path + "external_tools/blat -stepSize=5 -repMatch=2253 -minScore=5 -minIdentity=5 " + target_file + " " + assembly + " " + blat_output + "  > /dev/null 2>&1";
    // system(overlap_byBLAT.c_str());

    //blastn -subject  -query  -evalue 5  -word_size 4 -outfmt '7 qseqid qlen sseqid slen qcovs  qstart qend sstart send' -out 
    string overlap_BLAST = "blastn -evalue 5  -word_size 4 -outfmt '7 qseqid qlen sseqid slen qcovhsp  qstart qend sstart send' -subject " + target_file + " -query " + query_file + " -out " + blat_output + "  > /dev/null 2>&1";
    system(overlap_BLAST.c_str());

    string entry;
    vector<string> temp;

    ifstream blat_left(blat_output.c_str());

    while (!(getline(blat_left, entry).eof())) {

        if (entry[0] == '#') {
            continue;
        } else {
            temp.clear();
            utils.str_split(entry, temp, delimiter);
            //int query_start =
            //string strand = temp[8];
            int target_start = atoi(temp[7].c_str());
            int target_end = atoi(temp[8].c_str());
            //float score1 = atoi(temp[0].c_str());
            //  float final_score = score1 / query.size(); // score1 is number of matches

            if (flag == "left") {

                if (target_start < target.size() / 2) {

                    final = query + target.substr(target_end, target.size());
                }

            } else { // flag is right

                if (target_start > target.size() / 2) {

                    final = target.substr(0, target_start) + query;
                }
            }
            break;

        }
    }

    blat_left.close();
    remove(target_file.c_str());
    remove(query_file.c_str());
    remove(blat_output.c_str());

    return final;
}
/*string ROAST_extendContigs::overlap_mergeCAP3(string query, string target, string flag, int min_overlap_TH) {
    
    string overlap = "", target_region = "";
    string final = "";
    unsigned int found;
    boost::smatch present;
    int st_overlap_pos, j, en_overlap_pos;
    int overalapble_pattern = float(query.size() * 75 )/100;
    
    // const string pattern;
    bool isMatchFound, keep_checking = true;
    int a = 0, seed_size = 10, seed_limit = 5;
    // found indrx start from 1 and string index from 0
    // query is cap3 contig, target is trinity contig
    if (flag == "left") {
        
        while (keep_checking == true && seed_size >= seed_limit) {

            if (target.size() > overalapble_pattern)
                target_region = target.substr(0, overalapble_pattern); //target.size() / 2);
            else
                target_region = target.substr(0, target.size() / 2);
            
            string pattern = query.substr(query.size() - seed_size, query.size());

            boost::regex regexPattern(pattern, boost::regex::extended);
            isMatchFound = boost::regex_search(target_region, present, regexPattern);

            if (isMatchFound) {
                found = present.position(a);
                st_overlap_pos = found + seed_size;
                j = query.size() - (seed_size + 1);
                // cout << query[j] << endl;
                overlap = pattern;
                int i = found - 1;
                
                while (i >= query.size() - overalapble_pattern){ // 0) {
                    //     cout << target[i] << " " << query[j] << endl;
                    if (target[i] == query[j]) {
                        overlap = target[i] + overlap;
                        found = i;
                        j--;
                        i--;

                    } else break;
                }
                if (found <= 5 || overlap.size() >= min_overlap_TH) {
                    
                    std::transform(overlap.begin(), overlap.end(), overlap.begin(), ::tolower);
                    //cout << query.substr(0, j + 1) + overlap + target.substr(st_overlap_pos, target.size()) << endl;
                    final = query.substr(0, j + 1) + overlap + target.substr(st_overlap_pos, target.size());
                    keep_checking = false;
                    
                } else {
                    keep_checking = true;
                    seed_size--;
                }
            } else if (!isMatchFound) {
                keep_checking = true;
                seed_size--;
            }
        }

    }
    if (flag == "right") {
        
        while (keep_checking == true && seed_size >= seed_limit) {
            // found = target.find(query.substr(0, 5));
            string pattern = query.substr(0, seed_size);
            string target_region = target.substr(target.size() - overalapble_pattern, target.size() - 1);
            //  std::size_t found = target_region.find(query.substr(pattern));
            boost::regex regexPattern(pattern, boost::regex::extended);
            isMatchFound = boost::regex_search(target, present, regexPattern);
            
            if (isMatchFound) { // check 5 bases
                found = present.position(a);
                en_overlap_pos = found;
                j = seed_size;
                int i = found + seed_size;
                overlap = pattern;
                
                while (i < target.size() - overalapble_pattern) {

                    if (target[i] == query[j]) {
                        overlap = overlap + target[i];
                        found = i;
                        j++;
                        i++;
                    } else break;
                }
                if (found > target.size() - 5 || overlap.size() >= min_overlap_TH) {
                    
                    std::transform(overlap.begin(), overlap.end(), overlap.begin(), ::tolower);
                    //cout << target.substr(0, en_overlap_pos) + overlap + query.substr(j, query.size()) << endl;
                    final = target.substr(0, en_overlap_pos) + overlap + query.substr(j, query.size());
                    keep_checking = false;
                    
                } else {
                    keep_checking = true;
                    seed_size--;
                }
            } else if (!isMatchFound) {
                keep_checking = true;
                seed_size--;
            }
        }
    }
    return final;
}*/
int main(int argc, char** argv) {
    
 if (argc < 2) {
         cout << "No arguments found" << endl;
         exit(0);
     } 
      vector <string> args(argv, argv + argc);
      
      string bam_file = args[1];
      
      string in_file = args[2];
      
      string log = args[3];
      
      string num = args[4];
      
      string path_inter = args[5];
        
      string exe_path = args[6];

      int min_overlap_TH  = atoi(args[7].c_str()); //10;//
      
      int min_unmappead_reads_CAP3  = atoi(args[8].c_str()); //3

      int min_CAP3_ext = atoi(args[9].c_str());
      
      int read_length = atoi(args[10].c_str());

      string tool_name_flag = args[11];

      ROAST_extendContigs tt;
      tt.cap3_extension(bam_file, in_file, log, num, path_inter, exe_path, min_overlap_TH, min_unmappead_reads_CAP3, min_CAP3_ext, read_length, tool_name_flag);
      
     //cout <<"ROAST_extendContigs returned successfully " << endl;
      
       return 0;
}

