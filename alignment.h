/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   alignment.h
 * Author: madiha
 *
 * Created on April 26, 2019, 6:37 PM
 */

#ifndef ALIGNMENT_H
#define ALIGNMENT_H
#include<map>
#include <string>
#include <vector>
using namespace std;

class alignment {
public:
    alignment();
    alignment(const alignment& orig);
    virtual ~alignment();
    
    string align_reads(const string& reference_file, const string& left_reads_file, const string& right_reads_file);//, string cov_file) ;
    string align_reads(const string& reference_file, const string& left_reads_file, const string& right_reads_file, string cov_minimap, string cov_hisat);
    void create_indexes_minimap(const string& reference_file);
    void create_indexes_hisat(const string& reference_file);
    void calculate_coverage(string sam_sorted_filt, string cov_file);
    void remove_indexes(const string& reference_file);
    map<string, string>  extract_fasta_data(string fastafile);
    //void remove_indexes(const string& reference_file);
private:

};

#endif /* ALIGNMENT_H */

