/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ROAST_extendContigs.h
 * Author: madiha
 *
 * Created on December 23, 2020, 4:03 PM
 */

#ifndef ROAST_EXTENDCONTIGS_H
#define ROAST_EXTENDCONTIGS_H
#include<map>
#include <string>
#include <set>
#include <vector>
#include <list>
#include <fstream> 

using namespace std;

class ROAST_extendContigs {
public:
    ROAST_extendContigs();
    ROAST_extendContigs(const ROAST_extendContigs& orig);
    virtual ~ROAST_extendContigs();
    void cap3_extension(string bam_file, string in_file,  string &log, string num, string path_inter, string exe_path, int min_overlap_TH, int min_unmappead_reads_CAP3, int min_CAP3_ext, int read_length, string tool_name);
    string left_cap3_assembly(string bam_file, string debugFile,string unmapped_reads, string unmapped_mates, string fa_file, string assembly, string contig_name, string exe_path, int min_overlap_TH, int min_unmappead_reads_CAP3);
    string right_cap3_assembly(string bam_file, string debugFile, string unmapped_reads, string unmapped_mates, string fa_file, string assembly, string contig_name, string exe_path , string contig_length, int min_overlap_TH, int min_unmappead_reads_CAP3);
    string overlap_mergeCAP3(string query, string target, string flag, int min_overlap_TH, string path_inter, string exe_path, string num);
    
    
private:

};

#endif /* ROAST_EXTENDCONTIGS_H */

