/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   alignment.cpp
 * Author: madiha
 * 
 * Created on April 26, 2019, 6:37 PM
 */

#include <fstream>
#include "alignment.h"
#include "utils.h"
#include "iostream"
using namespace std;

alignment::alignment() {
}

alignment::alignment(const alignment& orig) {
}

alignment::~alignment() {
}

string alignment::align_reads(const string& reference_file, const string& left_reads_file, const string& right_reads_file) {//, string cov_file) {
    utils utils;

    stringstream ss;
    string thread_sam;

    if (allowed_threads > 1) {
        ss.str(string());
        ss << allowed_threads - 1;
        thread_sam = ss.str();
        ss.str(string());
    } else
        thread_sam = allowed_threads;


    string sam_file = reference_file.substr(0, reference_file.length() - 6) + ".sam";
    string bam_file = reference_file.substr(0, reference_file.length() - 6) + ".bam";
    string bam_sorted = reference_file.substr(0, reference_file.length() - 6) + ".sorted.bam"; //position sorted
    string bam_sorted_filt = bam_sorted.substr(0, bam_sorted.length() - 4) + ".filtered.bam";
    string bam_sorted_filt_temp = bam_sorted.substr(0, bam_sorted.length() - 4) + ".minimap.temp.bam";

    string command_minimap2, command_bam_sort, command_sam_to_bam, command_bam_to_cov, command_index;

    command_minimap2 = exe_path + "external_tools/minimap2 -Y -ax sr " + reference_file + " " + left_reads_file + " " + right_reads_file + " -o " + sam_file + " -t " + threads + "--secondary=no > /dev/null 2>&1";
    system(command_minimap2.c_str());

    command_sam_to_bam = "samtools view -bS " + sam_file + "  -o " + bam_file + " -@ " + thread_sam + "  > /dev/null 2>&1";
    system(command_sam_to_bam.c_str());

    command_bam_sort = "samtools sort " + bam_file + " -o " + bam_sorted + " -@ " + threadsForSamSortS + " -m " + MemForSamSortS + "M > /dev/null 2>&1";
    system(command_bam_sort.c_str());

    utils.filter_sam(bam_sorted, bam_sorted_filt_temp);
    
    // to filter supplementary alignments
    string command_filter = "samtools view -b -F 2048 " + bam_sorted_filt_temp + " -o " + bam_sorted_filt + " -@ " + thread_sam + "  > /dev/null 2>&1";
    system(command_filter.c_str());

    command_index = "samtools index " + bam_sorted_filt + " > /dev/null 2>&1";
    system(command_index.c_str());

    //remove indexes
   // remove_indexes(reference_file);
    utils.remove_file(sam_file);
    utils.remove_file(bam_file);
    utils.remove_file(bam_sorted);
    utils.remove_file(bam_sorted_filt_temp);

    //    return the sam object
    return bam_sorted_filt;
}

string alignment::align_reads(const string& reference_file, const string& left_reads_file, const string& right_reads_file, string cov_minimap, string cov_hisat) {
    utils utils;
    stringstream ss;
    ss << allowed_threads;
    threads = ss.str();
    string thread_sam;

    if (allowed_threads > 1) {
        ss.str(string());
        ss << allowed_threads - 1;
        thread_sam = ss.str();
        ss.str(string());
    } else
        thread_sam = allowed_threads;


    size_t find = reference_file.find_last_of(".");
    string sam_file = reference_file.substr(0, find) + ".sam";
    string bam_file = reference_file.substr(0, find) + ".bam";
    string bam_sorted = reference_file.substr(0, find) + ".sorted.bam"; //position sorted

    find = bam_sorted.find_last_of(".");
    string bam_sorted_minimap = bam_sorted.substr(0, find) + ".minimap.bam";
     string bam_sorted_minimap_temp = bam_sorted.substr(0, find) + ".minimap.temp.bam";
     
    string bam_sorted_hisat = bam_sorted.substr(0, find) + ".hisat2.bam";
    string merged_bam_file = bam_file.substr(0, find) + ".merged.bam";

    string bam = bam_file.substr(0, find) + ".bam";

    string command_minimap, command_bam_sort, command_sam_to_bam, command_bam_to_cov, command_cufflinks, command_mergeBam;

    string command_minimap2 = exe_path + "external_tools/minimap2  -Y -ax sr " + reference_file + " " + left_reads_file + " " + right_reads_file + " -o " + sam_file + " -t " + threads + " --secondary=no  > /dev/null 2>&1";
    system(command_minimap2.c_str());

    command_sam_to_bam = "samtools view -bS -F 4  " + sam_file + " -o " + bam + " -@ " + thread_sam + "  > /dev/null 2>&1";
    system(command_sam_to_bam.c_str());

    command_bam_sort = "samtools sort " + bam + " -o " + bam_sorted  + " -@ " + threadsForSamSortS + " -m " + MemForSamSortS + "M > /dev/null 2>&1";
    system(command_bam_sort.c_str());
    
    utils.filter_sam(bam_sorted, bam_sorted_minimap_temp);
    
     // to filter supplementary alignments
    string command_filter = "samtools view -b -F 2048 " + bam_sorted_minimap_temp + " -o " + bam_sorted_minimap + " -@ " + thread_sam + "  > /dev/null 2>&1";
    system(command_filter.c_str());
     
    string command_index = "samtools index " + bam_sorted_minimap + " > /dev/null 2>&1";
    system(command_index.c_str());

    command_bam_to_cov = "samtools depth -aa " + bam_sorted_minimap + " > " + cov_minimap + "  2> /dev/null";
    system(command_bam_to_cov.c_str());

    utils.remove_file(sam_file);
    utils.remove_file(bam);
    utils.remove_file(bam_sorted);
    utils.remove_file(bam_sorted_minimap_temp);


    //  hisat
    create_indexes_hisat(reference_file);
    string command_hisat2 = exe_path + "external_tools/hisat2 --downstream-transcriptome-assembly -f -x " + reference_file + " -q -1" + left_reads_file + " -2 " + right_reads_file + " -S " + sam_file + " --dta-cufflinks -p " + threads + "  > /dev/null 2>&1";
    system(command_hisat2.c_str());

    command_sam_to_bam = "samtools view -bS -F 4  " + sam_file + " -o " + bam + " -@ " + thread_sam + "  > /dev/null 2>&1";
    system(command_sam_to_bam.c_str());

    
    command_bam_sort = "samtools sort " + bam + " -o " + bam_sorted + " -@ " + threadsForSamSortS + " -m " + MemForSamSortS + "M > /dev/null 2>&1";
    system(command_bam_sort.c_str());

    utils.filter_sam(bam_sorted, bam_sorted_hisat);
    
    command_bam_to_cov = "samtools depth -aa " + bam_sorted_hisat + " > " + cov_hisat + "  2> /dev/null";
    system(command_bam_to_cov.c_str());


    // remove indexes
    remove_indexes(reference_file);
    utils.remove_file(sam_file);
    utils.remove_file(bam);
    utils.remove_file(bam_file);
    utils.remove_file(bam_sorted);

    //  cufflink intron-exon splice site
    size_t found = merged_bam_file.find_last_of("/\\");
    command_cufflinks = exe_path + "external_tools/cufflinks " + bam_sorted_hisat + " -o " + merged_bam_file.substr(0, found) + " --library-type fr-firststrand -p " + threads + " > /dev/null 2>&1";
    system(command_cufflinks.c_str());


    return bam_sorted_minimap;
}

void alignment::create_indexes_minimap(const string& reference_file) {


    stringstream ss;
    ss << allowed_threads;
    threads = ss.str();

    ss.str(string());
    ss << allowed_threads - 1;
    string thread_sam = ss.str();
    ss.str(string());
    string command = exe_path + "external_tools/minimap2-build -f " + reference_file + " " + reference_file + " -p " + threads + "  > /dev/null 2>&1"; //tmp_file;
    system(command.c_str());
}

void alignment::create_indexes_hisat(const string& reference_file) {
    //create the index
    string command = exe_path + "external_tools/hisat2-build -f " + reference_file + " " + reference_file + " -p " + threads + " > /dev/null 2>&1"; //tmp_file;
    system(command.c_str());
}

void alignment::remove_indexes(const string& reference_file) {
    utils utils;
    // remove index file
    utils.remove_file(reference_file + ".1.ht2");
    utils.remove_file(reference_file + ".2.ht2");
    utils.remove_file(reference_file + ".3.ht2");
    utils.remove_file(reference_file + ".4.ht2");
    utils.remove_file(reference_file + ".5.ht2");
    utils.remove_file(reference_file + ".6.ht2");
    utils.remove_file(reference_file + ".7.ht2");
    utils.remove_file(reference_file + ".8.ht2");

}

map<string, string> alignment::extract_fasta_data(string fastafile) {

    AllFasta_data.clear(); // empty before updating again
    string header = "", entry;
    stringstream seq;
    map<string, string> fasta_data;

    ifstream fasta_file;
    fasta_file.open(fastafile.c_str());

    while (!(getline(fasta_file, entry).eof())) {

        if (entry[0] == '>') { // contig numbers + 1
            fasta_data[header] = seq.str(); // global variable
            header = entry.substr(1, entry.size());
            seq.str(string());
        } else
            seq << entry;
    }
    fasta_data[header] = seq.str(); // global variable
    fasta_file.close();

    return fasta_data;
}