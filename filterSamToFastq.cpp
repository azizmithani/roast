/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   filterSamToFastq.cpp
 * Author: madiha
 * 
 * Created on April 26, 2019, 5:39 PM
 */

#include "filterSamToFastq.h"
#include "utils.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <set>

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAlignment.h"
#include "api/BamMultiReader.h"
#include "api/SamHeader.h"

using namespace BamTools;

using namespace std;

filterSamToFastq::filterSamToFastq() {
}

filterSamToFastq::filterSamToFastq(const filterSamToFastq& orig) {
}

filterSamToFastq::~filterSamToFastq() {
}

void filterSamToFastq::filterSam(string bam_in, string bam_out, string fastq_first_in, string fastq_first_out, string fastq_sec_in, string fastq_sec_out) {

    if (fastq_first_in.find("filt") != string::npos) {
        remove_prev_fastq(fastq_first_in, fastq_sec_in);
    }
    BamReader reader;
    bool pad;
  //  int FLAG_READ_UNMAPPED = 4;
   // int FLAG_MATE_UNMAPPED = 8;
    int count = 0;
    set <string> reads_to_keep;


    const SamHeader header = reader.GetHeader(); // returns header data object
    // Opens a bam file
    if (!reader.Open(bam_in)) {
        cerr << "Could not open BAM file." << endl;
        exit(0);
    }
    // RefData read_id, mate_id;
    vector<RefData> ref;
    BamTools::RefData ref_contig;
    ref = reader.GetReferenceData();
    BamWriter writer;

    if (!writer.Open(bam_out, header, ref)) {
        cerr << "Could not open output BAM file" << endl;
        exit(0);
        // return;
    }

    BamTools::BamAlignment al;
    while (reader.GetNextAlignment(al)) {
        string read_name = al.Name;
        if (al.RefID != -1 && al.MateRefID != -1) {
            //      cigar = al.CigarData;
            ref_contig = ref.at(al.RefID);
            int len_read_contig = ref_contig.RefLength;

            if (reads_to_keep.find(read_name) != reads_to_keep.end()) // read already added, add its mate too(mapped or unmapped both)
            {
                writer.SaveAlignment(al);
                //count++;
                reads_to_keep.insert(read_name);
            } else { //not already added? 
                //check both read and mate at the same time but add only read, mate will be added in first condition then
                if (!(al.IsMapped()) || (!(al.IsMateMapped()))) { // add all unmapped reads, mate will be added in first check then
                    writer.SaveAlignment(al);
                    reads_to_keep.insert(read_name);
                    //count++;

                } else { // if both are mapped and mate ref id is not -1 
                    BamTools::RefData mate_contig = ref.at(al.MateRefID);
                    int len_mate_contig = mate_contig.RefLength;

                    if ((al.Position <= left_edge_boundary) || (al.Position >= len_read_contig - right_edge_boundary) || (al.MatePosition <= left_edge_boundary) || (al.MatePosition >= len_mate_contig - right_edge_boundary)) {

                        writer.SaveAlignment(al);
                        reads_to_keep.insert(read_name);
                        //count++;
                    } else
                        continue;
                }
            }
        }
        else {
            writer.SaveAlignment(al);
            reads_to_keep.insert(read_name);
        }
            
    }
    reader.Close();
    writer.Close();
    remove(bam_in.c_str());
    filterFastq(bam_out, fastq_first_out, fastq_sec_out);
    //filterFastq(reads_to_keep, fastq_sec_in, fastq_sec_out);

    //cout<<count<<endl;

}

//void filterSamToFastq::filterFastq(set<string> reads_to_keep, string fastq_in, string fastq_out ){
void filterSamToFastq::filterFastq(string bam_out, string fastq_first_out, string fastq_sec_out ){

    //picard SamToFastq I= AT943.Filter4fastq.bam FASTQ=/AT943.Filter1.fastq SECOND_END_FASTQ= /AT943.Filter2.fastq VALIDATION_STRINGENCY=SILENT
    string command_picard_bam2fastq = "java -jar " + exe_path + "external_tools/picard.jar SamToFastq I=" + bam_out + " FASTQ=" + fastq_first_out + " SECOND_END_FASTQ=" + fastq_sec_out + " VALIDATION_STRINGENCY=SILENT  > /dev/null 2>&1" ;
    system(command_picard_bam2fastq.c_str());
    
}
void filterSamToFastq::remove_prev_fastq(const string& fastq1, const string& fastq2){
    utils utils;
    // remove index file
    utils.remove_file(fastq1);
    utils.remove_file(fastq2);
}