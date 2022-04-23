/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   filterSamToFastq.h
 * Author: madiha
 *
 * Created on April 26, 2019, 5:39 PM
 */
#include<map>
#include <string>
#include <vector>
#include <set>
using namespace std;
#ifndef FILTERSAMTOFASTQ_H
#define FILTERSAMTOFASTQ_H

class filterSamToFastq {
public:
    filterSamToFastq();
    filterSamToFastq(const filterSamToFastq& orig);
    virtual ~filterSamToFastq();
    void filterSam (string bam_in, string bam_out, string fastq_first_in, string fastq_first_out, string fastq_sec_in, string fastq_sec_out);
    void filterFastq(string bam_out,  string fastq_first_out, string fastq_sec_out);
    void remove_prev_fastq(const string& fastq1, const string& fastq2) ;
private:

};

#endif /* FILTERSAMTOFASTQ_H */

