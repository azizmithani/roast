 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   mis_assembly_chimera.cpp
 * Author: madiha
 * 
 * Created on December 29, 2019, 10:11 PM
 */

#include "mis_assembly_chimera.h"
#include "global.h"
#include "utils.h"
#include <sstream>
#include <iostream>
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
#include <list>
#include <math.h>   
#include <map>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAlignment.h>
#include <boost/regex.hpp>
#include <boost/any.hpp>
//#include <bamtools/api/BamAlignment.h>

using namespace BamTools;
using namespace boost;

mis_assembly_chimera::mis_assembly_chimera() {
}
mis_assembly_chimera::mis_assembly_chimera(const mis_assembly_chimera& orig) {
}

mis_assembly_chimera::~mis_assembly_chimera() {
}

int sc_no_th = 3;
int PROBLEM_MISSING_SEQ = 0;
int PROBLEM_MIS_ASSEMBLY = 1;
int PROBLEM_UNSUPPORTED_INSERTION = 2;
int PROBLEM_SELF_CHIMERA = 3;
int PROBLEM_MULTIGENE_CHIMERA = 4;
int PROBLEM_TRANSLOCATION = 5;
int PROBLEM_INVERSION = 6;
int PROBLEM_INVERSION_TRANSLOCATION = 7;
int PROBLEM_INVERSION_TRANSLOCATION_TRANSLOCATION = 8;
int PROBLEM_INVERSION_TRANSLOCATION_SELFCHIMERA = 9;
int PROBLEM_INVERSION_TRANSLOCATION_INVERSION_TRANSLOCATION = 10;
int PROBLEM_INVERSION_TRANSLOCATION_MULTIGENECHIMERA = 11;
int PROBLEM_INVERSION_TRANSLOCATIONandSELFCHIMERA = 12;

multimap <string, int>::iterator merged_contigs_itr;
int start_overlap_merged, end_overlap_merged;

// change TH to 100 for exempting exons

int mis_assembly_chimera::extract_chimera_mis_assembly(string bam_in, string coverage_file_hisat, string coverage_file_bowtie, string fastafile, string processed_fastafile, string cufflink_file , string mis_assembly_chimera_log) {
    /*extract coverage and splice-junction data*/
    utils utils;
    //SC_support_TH = SC_support_TH -  35; // to make it 40 for this module
    std::size_t found = fastafile.find_last_of("/\\");
    string PATH = fastafile.substr(0, found);
    string softclip_chimera = path_inter + "/mis_assembled_chimeric_positions";

    map<string, list<int> > exon_boundary =  utils.cufflink_data(cufflink_file);

   // int split_count = filterPos_processFasta(fastafile, processed_fastafile, softclip_chimera, mis_assembly_chimera_log, exon_boundary);


    map<string, vector<int> > cov_hisat = utils.coverage_data_updated(coverage_file_hisat);
    map<string, vector<int> > cov_bowtie = utils.coverage_data_updated(coverage_file_bowtie);

    contig_count = cov_bowtie.size();
    map<string, list<int> >::iterator ex_bo_it;

    vector<int> sc_pat_positions;
    sc_pat_positions.reserve(10);
    list<int> ex_bo;
    //cout << "done extracting coverage" << endl;


    int cov_th, cov_th_bow, dip = 0, steep = 0;
    vector<int> only_cov, only_cov_bow;


    float avg_cov, avg_cov_bow;

    vector<int> gradual_dip_chimera, gradual_dip_chimera2, coverage_dip_chimera; //single_SC_all,  single_SC, 
    gradual_dip_chimera.reserve(10); gradual_dip_chimera2.reserve(10);  coverage_dip_chimera.reserve(10); // single_SC_all.reserve(1000); single_SC.reserve(100);
    
    utils.remove_file(softclip_chimera);
    ofstream outfile;
    outfile.open(softclip_chimera.c_str(), fstream::app);

    BamReader reader;
    BamTools::BamAlignment al;
    vector<CigarOp> cigar;
    BamTools::RefData ref_contig;
    int deletion = 0, insertion = 0;
    vector<int> clipsize, read_pos, gen_pos;
    int prv_ref_id = 0, cur_ref_id, LSC_extended_size, RSC_extended_size, len_read_contig;
    bool pad;
    string contig_name, fasta_seq; //read_direction, mate_direction, 

    map<int, string> umappedMate_reads, dist_mapped_reads, same_dir_readpair, oppo_dir_readpair;
    multimap<int, string> right_SC_pos_seq, left_SC_pos_seq;

    // Opens a bam file
    if (!reader.Open(bam_in)) {
        cerr << "Could not open BAM file." << endl;
        exit(0);
    }

    vector<RefData> ref;
    ref = reader.GetReferenceData();
    //cout << "got bam data" << endl;

    while (reader.GetNextAlignment(al)) {

        if (al.RefID >= 0) { // to avoid unaligned reads which gives -1 RefID
            ref_contig = ref.at(al.RefID);
            cur_ref_id = al.RefID;
            if (prv_ref_id == cur_ref_id) { // deal one contig at a time
                //                 contig_name = ref_contig.RefName;
                cigar = al.CigarData;
                prv_ref_id = cur_ref_id;

                /* get softclip data using bamtools
                 * read_pos -> position of softclip on read
                 * gen_pos -> position of softclip on gene (take it)
                 * clipsize -> size of clipped bases, if size and read pos are same means clipping in the left or start of read
                 *             if size and read pos are diff means clipping at the end of read (right clipping)
                 */
                if (al.GetSoftClips(clipsize, read_pos, gen_pos, pad) == true) {

                    len_read_contig = ref_contig.RefLength;
                    contig_name = ref_contig.RefName;

                    ex_bo_it = exon_boundary.find(contig_name);
                    if (ex_bo_it != exon_boundary.end()) {
                        ex_bo = ex_bo_it->second;
                    }
                     
                    only_cov.reserve(len_read_contig + 1);
                    only_cov_bow.reserve(len_read_contig + 1);

                    for (vector<CigarOp>::iterator itt = cigar.begin(); itt < cigar.end(); itt++) {
                        char type = itt->Type;
                        if (type == 'D' || type == 'N')
                            deletion = deletion + itt->Length; //extract from read position in case deletion found
                    }
                    if (clipsize.size() == 1) { // if one sided clip , go for further checking
         
                        if (clipsize[0] == read_pos[0] && gen_pos[0] > min_SCs_MCprocess) { // don't take less than 25 from start or 25 of end SC as it shows partial of broken not chimera. 
         
                               for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 
                                   
                                    if (gen_pos[0] >= *boundary - 5 && (gen_pos[0] <= *boundary + 5)) {
                                        goto out_of_SC_check;
                                    }
                                }
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
                            if (!(sc_seq.empty()))
                                left_SC_pos_seq.insert(pair<int, string>(gen_pos[0], sc_seq)); // sc pos and seq
                            sc_seq = "";

                           // cout << al.QueryBases.substr(0, read_pos[0]) << endl;
                            //left_SC_pos_seq.insert(pair<int, string>(gen_pos[0], al.QueryBases.substr(0, read_pos[0]))); // sc pos and seq

                        } else if (clipsize[0] != read_pos[0] && gen_pos[0] < len_read_contig - min_SCs_MCprocess) {

                            for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 

                                if (gen_pos[0] >= *boundary - 5 && (gen_pos[0] <= *boundary + 5)) {
                                    goto out_of_SC_check;
                                }
                            }

                            int sc_start = 0, sc_end = 0;
                            string base_qualities = al.Qualities;
                            string sc_seq = "";
                            sc_start = read_pos[0] - deletion;

                            for (int sc_pos = read_pos[0] - deletion; sc_pos <= al.QueryBases.length() - 1; sc_pos++) {

                                int quality = int(base_qualities.at(sc_pos)) - 33;

                                if (quality >= base_quality) { // ascii code - 33 (for phred score)
                                    sc_end = sc_pos;
                                    sc_seq = al.QueryBases.substr(sc_start, sc_end);
                                } else
                                    break;
                            }
                            if (!(sc_seq.empty()))
                                right_SC_pos_seq.insert(pair<int, string>(gen_pos[0] + 1, sc_seq));
                            sc_seq = "";

                            // cout << al.QueryBases.substr(read_pos[0] - deletion, al.QueryBases.length() - 1);
                            //right_SC_pos_seq.insert(pair<int, string>(gen_pos[0] + 1, al.QueryBases.substr(read_pos[0] - deletion, al.QueryBases.length() - 1)));
                        }
                    } else {// if (clipsize.size() > 1) // if two sided clips, check if opposite side clip is < = 3 

                        if ((clipsize[0] == read_pos[0] && clipsize[1] <= max_SCs_oppositeSide_MCprocess) && gen_pos[0] > min_SCs_MCprocess) { //left SC //(clipsize[0] == read_pos[0]) && clipsize[0] >= min_sc_len  && clipsize[1] <= 3 //&& clipsize[0] >= min_sc_len 20dec2020

                            for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 

                                if (gen_pos[0] >= *boundary - 5 && (gen_pos[0] <= *boundary + 5)) {
                                    goto out_of_LSC_check;
                                }
                            }
                              int sc_start = 0, sc_end = 0;
                            string base_qualities = al.Qualities;
                            string sc_seq = "";
                            for (int sc_pos = 0; sc_pos <= read_pos[0]; sc_pos++) {
                                
                                 int quality = int(base_qualities.at(sc_pos)) - 33;
                                 
                                if (quality >= base_quality) { // ascii code - 33 (for phred score)
                                    sc_start = 0;
                                    sc_end = sc_pos;
                                    sc_seq = al.QueryBases.substr(0, sc_end);
                        }
                                else 
                                    break;
                                    
                            }
                            if (!(sc_seq.empty()))
                                left_SC_pos_seq.insert(pair<int, string>(gen_pos[0], sc_seq)); // sc pos and seq
                            sc_seq = "";
                            
                            //cout << al.QueryBases.substr(0, read_pos[0]) << endl;
                            //left_SC_pos_seq.insert(pair<int, string>(gen_pos[0], al.QueryBases.substr(0, read_pos[0])));
                        }
out_of_LSC_check:
                        ;

                        if ((clipsize[1] != read_pos[1] && clipsize[0] <= max_SCs_oppositeSide_MCprocess) && gen_pos[1] < len_read_contig - min_SCs_MCprocess) { //right SC clipsize[1] >= min_sc_len &&  25dec2020

                            for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 

                                if (gen_pos[1] >= *boundary - 5 && (gen_pos[1] <= *boundary + 5)) {
                                    goto out_of_SC_check;
                                }
                            }
                            int sc_start = 0, sc_end = 0;
                            string base_qualities = al.Qualities;
                            string sc_seq = "";
                            sc_start = read_pos[1] - deletion;

                            for (int sc_pos = read_pos[1] - deletion; sc_pos <= al.QueryBases.length() - 1; sc_pos++) {

                                int quality = int(base_qualities.at(sc_pos)) - 33;

                                if (quality >= base_quality) { // ascii code - 33 (for phred score)
                                    sc_end = sc_pos;
                                    sc_seq = al.QueryBases.substr(sc_start, sc_end);
                                } else
                                    break;
                            }
                            if (!(sc_seq.empty()))
                                right_SC_pos_seq.insert(pair<int, string>(gen_pos[1] + 1, sc_seq));
                            sc_seq = "";

                            //cout << al.QueryBases.substr(read_pos[1] - deletion, al.QueryBases.length() - 1);

                            // right_SC_pos_seq.insert(pair<int, string>(gen_pos[1] + 1, al.QueryBases.substr(read_pos[1] - deletion, al.QueryBases.length() - 1)));
                        }
                    }
out_of_SC_check:
                    ;

                    clipsize.clear();
                    read_pos.clear();
                    gen_pos.clear();
                    deletion = 0, insertion = 0;
                }
              //  RefData read_id = ref.at(al.RefID);
              //  RefData mate_id = ref.at(al.MateRefID);

                /*commented 17feb 2021 as we are not using this information any more*/
               /*if (al.IsReverseStrand())
                    read_direction = "R";
                else // if (!(al.AlignmentFlag & FLAG_READ_REVERSE_STRAND))
                    read_direction = "F";
                if (al.IsMateReverseStrand()) //.AlignmentFlag & FLAG_MATE_REVERSE_STRAND)
                    mate_direction = "R";
                else //if (!(al.AlignmentFlag & FLAG_MATE_REVERSE_STRAND))
                    mate_direction = "F";

                 if (!(al.IsMapped()) || (!(al.IsMateMapped()))) {
                    umappedMate_reads[al.Position] = read_direction;
                }// because = also appear for unmapped reads with * at cigar as it is merged file of two alignment we need double check
                else if ((al.IsMapped()) && (al.IsMateMapped()) && ((mate_id.RefName.compare(read_id.RefName) != 0) && )) { // -5 to avoid being in same gene group, but it would be subjective if we change contig id
                    dist_mapped_reads[al.Position] = read_direction;
                    //   cout<<al.Name << " " << read_id.RefName <<" " << mate_id.RefName <<endl;
                } else if ((!(al.AlignmentFlag & FLAG_MATE_UNMAPPED) && (!(al.AlignmentFlag & FLAG_READ_UNMAPPED))) && read_direction == mate_direction) { // 2 consecutive pos for read and mate
                    same_dir_readpair[al.Position] = read_direction;
                    same_dir_readpair[al.MatePosition] = mate_direction;
                }*/
                //                else if((!(al.AlignmentFlag&FLAG_MATE_UNMAPPED) && (!(al.AlignmentFlag&FLAG_READ_UNMAPPED))) && read_direction == mate_direction ){ //for R2F1 or F2R1
                //                
                //                }
                //   ref_contig = ref.at(al.RefID); //get name of specific contig id; 
            } else {
                // contig_name = ref_contig.RefName;
                //28 oct, collect start and end indexes of merged contigs to avoid chimera identification there 

                merged_contigs_itr = merged_contigs.find(contig_name);

                if (merged_contigs_itr != merged_contigs.end()) {
                    start_overlap_merged = merged_contigs_itr->second;
                    merged_contigs_itr++;
                    end_overlap_merged = merged_contigs_itr->second;
                }

                if (contig_name.find(ROASTcap3_left_tag) == string::npos && contig_name.find(ROASTcap3_right_tag) == string::npos) { // to avoid changes and splitting in cap3 generated assemblies
                    std::map<string, vector<int> >::iterator itr_hisat = cov_hisat.find(contig_name);
                    if (itr_hisat != cov_hisat.end()) { //extract coverage of that contig
                        only_cov = itr_hisat->second;

                    for (vector<int>::iterator a = only_cov.begin(); a != only_cov.end(); ++a){ // 31sMay2021 a++) {
                        avg_cov = avg_cov + *a;
                    }
                    avg_cov = avg_cov / only_cov.size();
                    cov_th = (avg_cov_th * avg_cov) / 100;
                }

                std::map<string, vector<int> >::iterator itr_bowtie = cov_bowtie.find(contig_name);
                unsigned currentMax = 0;
                
                if (itr_bowtie != cov_bowtie.end()) { //extract coverage of that contig
                    only_cov_bow = itr_bowtie->second;
                    
                    for (vector<int>::iterator b = only_cov_bow.begin(); b != only_cov_bow.end(); ++b) {
                        
                        if (*b > currentMax) {
                            
                            currentMax = *b;
                        }
                        avg_cov_bow = avg_cov_bow + *b;
                    }
                    
                    avg_cov_bow = avg_cov_bow / only_cov_bow.size();
                    cov_th_bow = (avg_cov_th * avg_cov_bow) / 100;
                }
                
                max_coverage = currentMax;
                

                if (only_cov.size() >= ignore_short_seq) {
                    //fasta_seq = utils.extract_fasta(contig_name, fastafile);
                      
                        map<string, string>::iterator it_fasta;
                        it_fasta = AllFasta_data.find(contig_name);

                        fasta_seq = it_fasta->second; //utils.extract_fasta(ID1, fastafile);
                        
                        if(coverage_drop > 0 && coverage_drop < 1 && only_cov.size() >= 600) //because wind_size is 200 each and read length from start and end
                            coverage_dip_chimera = chimera_cov1(contig_name, only_cov, cov_th); // chimera by sharp coverage change
                        
                        if(win_diff_TH > 0 && win_diff_TH < 1 && only_cov.size() >= 600) //because wind_size is 200 each and read length from start and end
                            gradual_dip_chimera = gradual_change(contig_name, only_cov, cov_th); // chimera by gradual coverage change

                        //old criteria to determine single sc pattern
                    //single_SC_all = single_SC_pattern(contig_name, left_SC_base_pos, left_SC_seq, right_SC_base_pos, right_SC_seq, only_cov_bow, cov_th_bow); // chimera by single sided SC

                    //*******************
                    //single_SC = filter_singleSC(single_SC_all, umappedMate_reads, dist_mapped_reads);
                    //                   **************
                    sc_pat_positions = sc_pattern(contig_name, fasta_seq, left_SC_pos_seq, right_SC_pos_seq, only_cov_bow, softclip_chimera, PATH, bam_in);
                

                    /*map<string, list<int> >::iterator ex_bo_it = exon_boundary.find(contig_name);
                    if (ex_bo_it != exon_boundary.end()) {
                        //ex_bo = ex_bo_it->second;
                    }*/
                    std::map<int, int>::iterator it1;
                    std::map<string, map<int, int> >::iterator it2;

                    if (!coverage_dip_chimera.empty()) {
                        //                
                        for (int it_cov_chi = 0; it_cov_chi < coverage_dip_chimera.size(); it_cov_chi++) {
                            if (!coverage_dip_chimera.empty() && !gradual_dip_chimera.empty()) {
                                for (int it_window = 0; it_window < gradual_dip_chimera.size(); it_window++) {
                                    if (abs(coverage_dip_chimera[it_cov_chi] - gradual_dip_chimera[it_window]) <= consecutive_missAssembled_pos_dist) {
                                        if (only_cov[gradual_dip_chimera[it_window] ] >= only_cov[coverage_dip_chimera[it_cov_chi]]) {
                                            gradual_dip_chimera.erase(gradual_dip_chimera.begin() + it_window);
                                            it_window = it_window - 1;
                                        } else {
                                            coverage_dip_chimera.erase(coverage_dip_chimera.begin() + it_cov_chi);
                                            it_cov_chi = it_cov_chi - 1;
                                            goto loop_cov;
                                        }
                                    }
                                }
                            }

                            if (ex_bo_it != exon_boundary.end()) { // check splice-junction
                                for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) {
                                    if (coverage_dip_chimera[it_cov_chi] >= *boundary - 5 && (coverage_dip_chimera[it_cov_chi] <= *boundary + 5)) {
                                        coverage_dip_chimera.erase(coverage_dip_chimera.begin() + it_cov_chi);
                                        it_cov_chi = it_cov_chi - 1;
                                        steep++;
                                        goto loop_cov;
                                    }
                                }
                            }

                            if (!sc_pat_positions.empty() && !coverage_dip_chimera.empty()) { // check miss-assembled regions
                                for (vector<int>::iterator sc = sc_pat_positions.begin(); sc != sc_pat_positions.end(); sc++) {
                                    if (abs(coverage_dip_chimera[it_cov_chi] - *sc) <= consecutive_missAssembled_pos_dist) { // 100 because of read length to dilute effect of cov change due to SCs 
                                        coverage_dip_chimera.erase(coverage_dip_chimera.begin() + it_cov_chi);
                                        it_cov_chi = it_cov_chi - 1;
                                        steep++;
                                        goto loop_cov;
                                    }
                                }
                            }
                                if (!coverage_dip_chimera.empty()) {
                                    string type = chimera_type(coverage_dip_chimera[it_cov_chi], fasta_seq);

                                    //28oct21 avoid chimeras at the start or end positions of overlapped merged positions found in RI and merged by SCs
                                    if (merged_contigs_itr != merged_contigs.end()) {
                                        if (coverage_dip_chimera[it_cov_chi] > end_overlap_merged || coverage_dip_chimera[it_cov_chi] < start_overlap_merged) {
                                            //if (type.find("multigeneChimera") != string::npos) { // don't make multigeneChimera of merged contigs by RI and SCs

                                            outfile << contig_name << "\t" << type << "\t" << "sharp coverage change" << endl;
                                            //   }
                                            // } else {
                                            //  outfile << contig_name << "\t" << type << "\t" << "sharp coverage change" << endl;
                                            //  }
                                            // }
                                            //cout << "em in else" << endl;
                                        }
                                    }
                                    else
                                        outfile << contig_name << "\t" << type << "\t" << "sharp coverage change" << endl;
                                }
loop_cov:
                                    ;
                               // }
                        if (!gradual_dip_chimera.empty()) {

                                    gradual_dip_chimera2 = compare_gradual_cov_change(gradual_dip_chimera, ex_bo, sc_pat_positions); //with soft clips,keep intron check

                                    if (!gradual_dip_chimera2.empty()) {

                                        for (int it_window = 0; it_window < gradual_dip_chimera2.size(); it_window++) {

                                            string type = chimera_type(gradual_dip_chimera2[it_window], fasta_seq);

                                            //28oct21 avoid chimeras at the start or end positions of overlapped merged positions found in RI and merged by SCs
                                            if (merged_contigs_itr != merged_contigs.end()) {
                                                if (gradual_dip_chimera2[it_window] > end_overlap_merged || gradual_dip_chimera2[it_window] < start_overlap_merged) {
                                                    outfile << contig_name << "\t" << type << "\t" << "gradual coverage change" << endl;
                                                }
                                            } else
                                                outfile << contig_name << "\t" << type << "\t" << "gradual coverage change" << endl;
                                        }
                                    }
                                    gradual_dip_chimera2.clear();
                                }
                            }//end of first condition

                    }
                    else if (!gradual_dip_chimera.empty()) {

                                gradual_dip_chimera2 = compare_gradual_cov_change(gradual_dip_chimera, ex_bo, sc_pat_positions);

                                if (!gradual_dip_chimera2.empty()) {

                                    for (int it_window = 0; it_window < gradual_dip_chimera2.size(); it_window++) {

                                        string type = chimera_type(gradual_dip_chimera2[it_window], fasta_seq);

                                    //28oct21 avoid chimeras at the start or end positions of overlapped merged positions found in RI and merged by SCs
                                    if (merged_contigs_itr != merged_contigs.end()) {
                                        if (gradual_dip_chimera2[it_window] > end_overlap_merged || gradual_dip_chimera2[it_window] < start_overlap_merged) {

                                            outfile << contig_name << "\t" << type << "\t" << "gradual coverage change" << endl;
                                        }
                                    } else
                                        outfile << contig_name << "\t" << type << "\t" << "gradual coverage change" << endl;
                                        //outfile << contig_name << "\t" << gradual_dip_chimera2[it_window] << "\t" << "gradual coverage" << endl;
                                    }
                                }
                                gradual_dip_chimera2.clear();
                            }//end of second condition

                        } // ignore_short_seq

                } //to avoid splitting cap3 assemblies
                //   left_SC_base_pos.clear();
                dist_mapped_reads.clear();
                // single_SC_all.clear();
                umappedMate_reads.clear();
                //left_SC_seq.clear();
                // right_SC_base_pos.clear();
                //right_SC_seq.clear();
                //  single_SC.clear();
                ex_bo.clear();
                only_cov.clear();
                only_cov_bow.clear();
                coverage_dip_chimera.clear();
                gradual_dip_chimera.clear();
                sc_pat_positions.clear();
                right_SC_pos_seq.clear();
                left_SC_pos_seq.clear();
                prv_ref_id = cur_ref_id; // got new contig

                contig_name = ref_contig.RefName;
            }
        }
    }
        
    reader.Close();
    outfile.close();
    //f1.longest_ORF(softclip_chimera, fastafile, pot_File, overlap, cov_hisat, cov_bowtie);
    int split_count = filterPos_processFasta(fastafile, processed_fastafile, softclip_chimera, mis_assembly_chimera_log, exon_boundary);
    return split_count;
}

// update 11/11/19: old criteria to determine single sc

    string mis_assembly_chimera::chimera_type(int pos, string contig_seq) {
    utils utils;
    string sec_frag, query_frag, blast_found, blast_hit, chimera_type;
    int start, end;
    vector<string> temp;
    stringstream entry_new;
    string first_frag = contig_seq.substr(0, pos);
    if (pos + 1 < contig_seq.size())
        sec_frag = contig_seq.substr(pos + 1, contig_seq.size()); // fragment after
    else
        sec_frag = "";

    if (first_frag.size() < sec_frag.size()) { // take smallest fragment of contig to check for chimera or local mis-assembly
        query_frag = first_frag;
        first_frag = "";
        start = 0;
        end = pos;
    } else { // sec fragment is smaller
        query_frag = sec_frag;
        sec_frag = "";
        
        start = pos;
        end = contig_seq.length() -1;
    }
    blast_hit = perform_blast(query_frag, first_frag, sec_frag);
 /* bug removed: 10june2021*/
    temp.clear();
    utils.str_split(blast_hit, temp, delimiter);
    blast_found = temp[0];

    if (blast_found != "*") {
        chimera_type = "Self";
        if (blast_found == "+") {
            entry_new << "selfChimera" << "\t" << start << "-" << end << "\t" << "Forward" << "\t" << "Discard";
        } else {
            query_frag = utils.Rcomplement(query_frag);
            entry_new << "selfChimera" << "\t" << start << "-" << end << "\t" << "Reverse complement" << "\t" << "Discard";
        }
    } else {
        if (abs(end - start) >= ignore_short_seq) { // no cov check because its either one side of contig, not in the middle 
            entry_new << "multigeneChimera" << "\t" << start << "-" << end << "\t" << "Forward" << "\t" << "New contig";

        } else
            entry_new << "mis-assembly" << "\t" << start << "-" << end << "\t" << "Forward" << "\t" << "Discard";
    }
    chimera_type = entry_new.str();
    entry_new.str(string());
    return chimera_type;

}

vector<int> mis_assembly_chimera::single_SC_pattern(string contig_name, vector<int>left_SC_base_pos, vector<string>left_SC_seq, vector<int>right_SC_base_pos, vector<string>right_SC_seq, vector<int> cov, int cov_th) {
   utils utils;
    //int SC_support_TH = 75;
    //cout << "entering single-sided sc" << endl;
    std::map<int, int>::iterator itr1_cov;
    int last_count = 0, current_count = 0;
    int last_leftSC = 0, last_rightSC = 0;
    bool rem = false;
    float consective_leftSC = 0, consective_rightSC = 0;
    float previous_last_leftSC = 0, previous_last_rightSC = 0;
    
    vector<string> left_seq, right_seq;
    left_seq.reserve(max_coverage); right_seq.reserve(max_coverage);
    
    vector<int> single_rSC, single_lSC, single_SC;
    single_rSC.reserve(max_coverage); single_lSC.reserve(max_coverage); single_SC.reserve(max_coverage);
    
    for (int m = 1; m < left_SC_base_pos.size(); m++) {
        if (left_SC_base_pos[m] - left_SC_base_pos[m - 1] >= 5) //if at a distance
        {
            left_SC_base_pos[m - 1] = 0;
            left_SC_seq[m - 1] = "";
        } else {
            continue;
        }
    }
    for (int k = 1; k < right_SC_base_pos.size(); k++) {
        if (right_SC_base_pos[k] - right_SC_base_pos[k - 1 ] >= 5) {
            right_SC_base_pos[k - 1] = 0;
            right_SC_seq[k - 1] = " ";
        } else {
            continue;
        }
    }
    for (int mm = 1; mm < left_SC_base_pos.size(); mm++) { // make clusters of max blocks, filter by coverage and consecutive read patterns
        if (left_SC_base_pos[mm] != 0 && mm <= left_SC_base_pos.size() - 2) {
            if (left_SC_base_pos [mm + 1] - left_SC_base_pos[mm] == 0) // to get the position of max support
                last_leftSC = left_SC_base_pos[mm];
            consective_leftSC += 1;
            //     left_seq.push_back(left_SC_seq[mm]);
        } else {
            if (consective_leftSC > sc_no_th) {// && max_consec_leftSC < consective_leftSC) {
                //    if (previous_last_leftSC != last_leftSC) {
                previous_last_leftSC = last_leftSC;
                for (vector<int>::iterator c = cov.begin(); c != cov.end(); c++) {
                    int i = std::distance(cov.begin(), c) - 1; // 1 because cov index from 0
                    if (i == previous_last_leftSC && float(consective_leftSC / *c)*100 >= SC_support_TH) { // coverage check
                        //   cout << previous_last_leftSC << "  " << itr1_cov->first << " " << itr1_cov->second << endl;
                        if (single_lSC.size() > 0) {
                            current_count = consective_leftSC; // itr1_cov->second;
                            if (previous_last_leftSC - single_lSC.back() < 200) { // two consecutive sc patterns with distance
                                if (current_count > last_count) {
                                    single_lSC.pop_back();
                                    single_lSC.push_back(previous_last_leftSC);
                                    last_count = current_count;
                                }
                                //else do nothing
                            } else {
                                single_lSC.push_back(previous_last_leftSC);
                                last_count = current_count;
                            }

                        } else {
                            single_lSC.push_back(previous_last_leftSC);
                            last_count = consective_leftSC; //itr1_cov->second;
                        }
                        break;
                    }

                }
                //    last_cov=0, current_cov=0;
                //  }
            }
            //     left_seq.clear();
            last_leftSC = 0;
            consective_leftSC = 0;
        }
    }

    // extract single right soft clips
    for (int kk = 1; kk < right_SC_base_pos.size(); kk++) {
        //   cout<<right_SC_base_pos[kk]<<endl;
        if (right_SC_base_pos[kk] != 0 && kk <= right_SC_base_pos.size() - 2) {
            if (right_SC_base_pos [kk + 1] - right_SC_base_pos[kk] == 0) // to get the position of max support
                last_rightSC = right_SC_base_pos[kk];
            consective_rightSC += 1;
            // right_seq.push_back(right_SC_seq[kk]);
        } else {
            if (consective_rightSC > sc_no_th)// && max_consec_rightSC < consective_rightSC) //max cluster supported by  >4 SC reads
            {
                previous_last_rightSC = last_rightSC;
                for (vector<int>::iterator c = cov.begin(); c != cov.end(); c++) {
                    int i = std::distance(cov.begin(), c) - 1; // 1 because cov index from 0
                    if (i == previous_last_rightSC && float(consective_rightSC / *c)*100 >= SC_support_TH) {
                        //   cout << previous_last_rightSC << "  " << itr1_cov->first << " " << itr1_cov->second << endl;
                        if (single_rSC.size() > 0) {
                            current_count = consective_rightSC; //itr1_cov->second;
                            if (previous_last_rightSC - single_rSC.back() < 200) {
                                if (current_count > last_count) {
                                    single_rSC.pop_back();
                                    single_rSC.push_back(previous_last_rightSC);
                                    last_count = current_count;
                                }
                                //else do nothing
                            } else {
                                single_rSC.push_back(previous_last_rightSC);
                                last_count = current_count;
                            }

                        } else {
                            single_rSC.push_back(previous_last_rightSC);
                            last_count = consective_rightSC;
                        }
                        break;

                    }
                }
            }
            last_rightSC = 0;
            consective_rightSC = 0;
        }
    }

    if (single_lSC.size() > 0 && single_rSC.size() > 0) { //for 1 left SC cluster check all right soft clip clusters
        for (int a = 0; a < single_lSC.size(); a++) {
            for (int b = 0; b < single_rSC.size(); b++) {
                if (single_lSC[a] - single_rSC[b] <= 200 && single_lSC[a] - single_rSC[b] >= -200)//checked both ELP and SCchimera pattern
                {
                    single_rSC.erase(single_rSC.begin() + b);
                    b = b - 1;
                    rem = true;
                }
            }
            if (rem) {
                single_lSC.erase(single_lSC.begin() + a);
                a = a - 1;
                rem = false;
            }
            if (single_lSC.size() == 0)
                break;

        }
        for (int a = 0; a < single_lSC.size(); a++) {
            if (single_lSC[a] > st_end_th && single_lSC[a] < cov.size() - st_end_th && cov[single_lSC[a] - 1] <= cov_th) {
                single_SC.push_back(single_lSC[a]);
                //cout << contig_name << " " << single_lSC[a] << endl;
            }
        }
        for (int b = 0; b < single_rSC.size(); b++) {
            if (single_rSC[b] > st_end_th && single_rSC[b] < cov.size() - st_end_th && cov[single_rSC[b] - 1] <= cov_th) {
                single_SC.push_back(single_rSC[b]);
                // cout << contig_name << " " << single_rSC[b] << endl;
            }
        }
    } else {
        for (int a = 0; a < single_lSC.size(); a++) {
            if (single_lSC[a] > 100 && single_lSC[a] < cov.size() - st_end_th && cov[single_lSC[a] - 1] <= cov_th) {
                single_SC.push_back(single_lSC[a]);
                // cout << contig_name << " " << single_lSC[a] << endl;
            }
        }
        for (int b = 0; b < single_rSC.size(); b++) {
            if (single_rSC[b] > st_end_th && single_rSC[b] < cov.size() - st_end_th && cov[single_rSC[b] - 1] <= cov_th) {
                single_SC.push_back(single_rSC[b]);
                // cout << contig_name << " " << single_rSC[b] << endl;
            }
        }
    }
    return single_SC;

}

vector< int> mis_assembly_chimera::chimera_cov1(string contig_name, vector<int> coverage, int cov_th) //for identifying all drops except indels and exons
{
    utils utils;
    //cout << "entering sharp cov" << endl;
    vector<int> coverage_chimera, coverage_chimera_final;
    coverage_chimera.reserve(50); coverage_chimera_final.reserve(50);
    
    int window_TH = 100, size = 0, size1 = 0, setof = 50;
    
    vector<int> base_pos, ind;
    std::map<int, int>::iterator itr_cov;

    for (int a = 0; a < coverage.size(); a++) { //31st May 2021
        base_pos.push_back(a);
    }
    for (int b = setof; b < coverage.size() - setof; b++) { // 51 because coverage index pos from 
        if ((coverage[b] >= 6) && (coverage[b] > coverage[b - 1]) && (float(float(coverage[b - 1]) / float(coverage[b])) < coverage_drop)) {//sharp rise // 6 because we are not taking max cov <=5 
            ind.push_back(b - 1); // 
        } else if ((coverage[b - 1] >= 6) && (coverage[b - 1] > coverage[b]) && (float(float(coverage[b ]) / float(coverage[b - 1])) < coverage_drop)) //sharp fall
        {
            ind.push_back(b);
        }
    }
    if (ind.size() < 1) {
        //  cout << "nothing found" << endl;
    } else if (ind.size() == 1) {
        coverage_chimera.push_back(ind[0]);
    } else if (ind.size() == 2) {
        if ((ind[1] - ind[0]) <= indel_TH) {
            //do nothing
        } else if ((ind[1] - ind[0]) <= window_TH) {
            coverage_chimera.push_back(ind[0]);

        } else {
            coverage_chimera.push_back(ind[0]);
            coverage_chimera.push_back(ind[1]);
        }
    } else {
        if ((ind[1] - ind[0]) <= indel_TH) {
            //do nothing
            size1++;
        } else if ((ind[1] - ind[0]) <= window_TH) {
            coverage_chimera.push_back(ind[0]);
            size1++;
        } else {
            coverage_chimera.push_back(ind[0]);
            coverage_chimera.push_back(ind[1]);
        }
        for (int it = 2; it < ind.size() - 2; it++) {
            if (((ind[it] - ind[it - 1]) <= indel_TH) && ((ind[it + 1] - ind[it]) <= indel_TH)) { // check indel
                //    do nothing as its indel
                // coverage_chimera.push_back(base_pos[ind[it]]);
            } else if (((ind[it] - ind[it - 1]) <= window_TH) || ((ind[it + 1] - ind[it]) <= window_TH)) {
                if (((ind[it] - ind[it - 1]) <= window_TH) && (ind[it + 1] - ind[it]) > window_TH) {
                    if (coverage_chimera.size() > 0 && ind[it] - coverage_chimera[coverage_chimera.size() - 1] > window_TH)
                        coverage_chimera.push_back(ind[it]);
                    else if (coverage_chimera.size() == 0) coverage_chimera.push_back(ind[it]);
                    //    continue;
                } else {
                    if (coverage_chimera.size() > 0 && ind[it] - coverage_chimera[coverage_chimera.size() - 1] > window_TH)
                        coverage_chimera.push_back(ind[it]);
                    else if
                        (coverage_chimera.size() == 0) coverage_chimera.push_back(ind[it]);
                    it++;
                }
            } else {
                coverage_chimera.push_back(ind[it]);
            }
            size++;
        }
        if (size > 0 && size1 > 0) {
            if ((ind[ind.size() - 1] - ind[ind.size() - 2]) <= indel_TH) {
                //coverage_chimera.push_back(base_pos[ind[0]]);
                //do nothing
            } else if ((ind[ind.size() - 1] - ind[ind.size() - 2]) <= window_TH) {
                if (coverage_chimera.size() > 0 && ind[ind.size() - 1] - coverage_chimera[coverage_chimera.size() - 1] > window_TH)
                    coverage_chimera.push_back(ind[ind.size() - 1]);
                else if (coverage_chimera.size() == 0) coverage_chimera.push_back(ind[ind.size() - 1]);
            } else {
                if (coverage_chimera.size() > 0 && ind[ind.size() - 1] - coverage_chimera[coverage_chimera.size() - 1] > window_TH) {
                    coverage_chimera.push_back(ind[ind.size() - 2]);
                    coverage_chimera.push_back(ind[ind.size() - 1]);
                } else if (coverage_chimera.size() == 0) coverage_chimera.push_back(ind[ind.size() - 1]);
            }


            size = 0, size1 = 0;
        } else if (size == 0 && size1 > 0) { // to cater 333, 357, 478 type scenerios


            if (coverage_chimera.size() > 0 && ind[ind.size() - 1] - coverage_chimera[coverage_chimera.size() - 1] > window_TH)
                coverage_chimera.push_back(ind[ind.size() - 1]);
            else if (coverage_chimera.size() == 0) coverage_chimera.push_back(ind[ind.size() - 1]);
        }
    }
    // cout<< contig_name<<" in coverage2 "<<cov.size()<<" "<<endl;
    for (int it = 0; it < coverage_chimera.size(); it++) {
        if (coverage_chimera[it] >= st_end_th && coverage_chimera[it] <= coverage.size() - st_end_th && coverage[coverage_chimera[it]] <= cov_th) {
            coverage_chimera_final.push_back(coverage_chimera[it]);
        }
    }
    return coverage_chimera_final;
}

vector<int> mis_assembly_chimera::gradual_change(string contig_name, vector<int> coverage, int cov_th) {
    utils utils;
    //cout << "entering gradual" << endl;
    float left_count = 0, right_count = 0;
    float ratio;
    vector <string> temp;
    vector<float> pos;
    pos.reserve(100);
    
    vector<int> e_pos, g_pos, f_pos;
    e_pos.reserve(50); g_pos.reserve(50); f_pos.reserve(50);
    
    float diff_at_pos = 0.0, prev_diff_at_pos = 0.0;
    float pos_cov_left, pos_cov_right, pos_cov, window_TH = 100, fluc_TH = 20; //200
    bool second_round = false;
    int position = 0, prev_pos = 0, next_pos = 0, found_pos = 0, fluctuation = 0, size = 0, size1 = 0, read_length2 = read_length*3;
    int main_pos, p_pos, n_pos; //150, 299
    
    std::vector<int>::iterator itr_cov;
    //    for (itr_cov = coverage.begin(); itr_cov != coverage.end(); itr_cov++) {
    //        cov.push_back(itr_cov->second); // cout<<itr1->second<<endl;
    //    }

    for (int i = read_length; i < coverage.size() - read_length2; i++) { // coverage starts from 0 index pos instead 1
        int window1 = 0, window2 = 0;
        float avg1, avg2;
        int window3 = 0, window4 = 0;
        float avg3, avg4;
        for (int j = i; j < i + win_size; j++) //window size
        {
            window1 = window1 + coverage[j];
            window2 = window2 + coverage[j + win_size + 1]; // gap between window

        }
        avg1 = window1 / win_size;
        avg2 = window2 / win_size;
        next_pos = (i + win_size); // + gap);
        if (avg1 < avg2) {
            diff_at_pos = avg1 / avg2;
        } else {
            diff_at_pos = avg2 / avg1;
        }
        if (diff_at_pos < win_diff_TH && found_pos == 0) {
            position = next_pos;
            main_pos = next_pos;
            p_pos = main_pos - 1;
            n_pos = main_pos + 1;
            prev_pos = main_pos - 1;
            prev_diff_at_pos = diff_at_pos;
            found_pos++;
        } else if (found_pos > 0) { // roll down thing

            if (next_pos >= coverage.size() - (read_length + 50)) { //roll down select initial 150
                position = coverage.size() - (read_length + 52); //main_pos; 152
                pos_cov = coverage[prev_pos - 1];
                pos.push_back(diff_at_pos);
                g_pos.push_back(position);
                break;
            } else if (prev_pos <= (read_length + 50)) { //150
                position = read_length + 51; // main_pos; //151
                pos_cov = coverage[position - 1];
                pos.push_back(diff_at_pos);
                g_pos.push_back(position);
                i = i + win_size;
                found_pos = 0;
                continue;

            } else if ((avg1 > avg2) && (coverage[next_pos - 1] <= coverage[position - 1])) {
                position = next_pos;
                //prev_diff_at_pos = diff_at_pos;
                fluctuation = 0;
                second_round = true;
            } else if ((avg1 > avg2) && (coverage[next_pos - 1] - coverage[position - 1] <= 1 && coverage[next_pos - 1] - coverage[position - 1] > 0) && fluctuation < fluc_TH) {
                fluctuation++;
            } else if ((avg2 > avg1) &&(coverage[prev_pos - 1] <= coverage[position - 1])) {
                position = prev_pos;
                prev_pos--;
                fluctuation = 0;
                second_round = true;
                //       cout<<cov[position - 1] << " "<<cov[prev_pos - 1]<<endl; 
            } else if ((avg2 > avg1) && (coverage[prev_pos - 1] - coverage[position - 1] <= 1 && coverage[prev_pos - 1] - coverage[position - 1] > 0) && fluctuation < fluc_TH) {
                fluctuation++;
                prev_pos--;
            } else if (avg1 > avg2 && (coverage[next_pos - 1] - coverage[position - 1] >= 2) && fluctuation < 10) {
                fluctuation++;
            } else if (avg2 > avg1 && (coverage[prev_pos - 1] - coverage[position - 1] >= 2) && fluctuation < 10) {
                fluctuation++;
                prev_pos--;
            } else {
                i = i + win_size;
                found_pos = 0;
                fluctuation = 0;

                for (int k = position - 50; k <= position + 50; k++) {
                    if (position > 200 && position < coverage.size() - 200 && coverage[k] <= cov_th && coverage[k] < coverage[position - 1])
                        position = k;
                }
                pos_cov = coverage[position - 1];
                pos.push_back(diff_at_pos);
                g_pos.push_back(position);
            }

        }
    }
    if (pos.size() < 1) {
    } else if (g_pos.size() == 1) {
        //  cout << "nothing found" << endl;
        e_pos.push_back(g_pos[0]);
    } else if (g_pos.size() == 2) {
        if ((g_pos[1] - g_pos[0]) <= window_TH) {
            if (coverage[g_pos[0] - 1] < coverage[g_pos[1] - 1]) e_pos.push_back(g_pos[0]);
            else e_pos.push_back(g_pos[1]);
        } else {
            e_pos.push_back(g_pos[0]);
            e_pos.push_back(g_pos[1]);
        }
    } else {
        for (int it = 1; it < pos.size() - 1; it++) {
            if (((g_pos[it] - g_pos[it - 1]) <= window_TH) || ((g_pos[it + 1] - g_pos[it]) <= window_TH)) { // check indel
                //    cout << "indel" << endl;
                //                pos.erase(pos.begin()+ it);
                //                it=it-1;
                if (((g_pos[it] - g_pos[it - 1]) <= window_TH) && (g_pos[it + 1] - g_pos[it]) > window_TH) {
                    // g_pos.push_back(pos[it]);
                    // cout << cov[g_pos[it] - 1] << " " << cov[g_pos[it - 1] - 1] << endl;
                    if (coverage[g_pos[it] - 1] < coverage[g_pos[it - 1] - 1]) {
                        if (!(e_pos.empty())) {
                            if (g_pos[it] - e_pos[e_pos.size() - 1] > window_TH) {
                                e_pos.push_back(g_pos[it]);
                            }
                        } else {
                            e_pos.push_back(g_pos[it]);
                        }
                    } else {
                        if (!(e_pos.empty())) {
                            if (g_pos[it - 1] - e_pos[e_pos.size() - 1] > window_TH) {
                                e_pos.push_back(g_pos[it - 1]);
                            }
                        } else {
                            e_pos.push_back(g_pos[it - 1]);
                        }

                    }

                    // continue;
                } else if (((g_pos[it] - g_pos[it - 1]) > window_TH) && (g_pos[it + 1] - g_pos[it]) <= window_TH) {
                    if (!(e_pos.empty())) {
                        if (g_pos[it - 1] - e_pos[e_pos.size() - 1] > window_TH) {
                            e_pos.push_back(g_pos[it - 1]);
                        }
                    } else {
                        e_pos.push_back(g_pos[it - 1]);
                    }
                    if (coverage[g_pos[it] - 1] < coverage[g_pos[it + 1] - 1]) {
                        if (g_pos[it] - e_pos[e_pos.size() - 1] > window_TH) {
                            //pos.erase(pos.begin()-1);
                            e_pos.push_back(g_pos[it]);
                        }
                    } else {
                        if (g_pos[it + 1] - e_pos[e_pos.size() - 1] > window_TH) {
                            //pos.erase(pos.begin()-1);
                            e_pos.push_back(g_pos[it + 1]);
                        }
                    }
                    it++;
                } else {


                    if (coverage[g_pos[it] - 1] < coverage[g_pos[it + 1] - 1] && coverage[g_pos[it] - 1] < coverage[g_pos[it - 1] - 1]) { //check current
                        if (!(e_pos.empty())) {
                            if (g_pos[it] - e_pos[e_pos.size() - 1] > window_TH) {
                                //pos.erase(pos.begin()-1);
                                e_pos.push_back(g_pos[it]);
                            }
                        } else {
                            e_pos.push_back(g_pos[it]);
                        }
                    } else if (coverage[g_pos[it - 1] - 1] < coverage[g_pos[it] - 1] && coverage[g_pos[it - 1] - 1] < coverage[g_pos[it + 1] - 1]) { //check previous
                        if (!(e_pos.empty())) {
                            if (g_pos[it - 1] - e_pos[e_pos.size() - 1] > window_TH) {
                                //pos.erase(pos.begin()-1);
                                e_pos.push_back(g_pos[it - 1]);
                            }
                        } else {
                            e_pos.push_back(g_pos[it - 1]);
                        }
                    } else { // check next
                        if (g_pos[it + 1] - g_pos[it - 1] > window_TH) {
                            if (!(e_pos.empty())) {
                                if (g_pos[it - 1] - e_pos[e_pos.size() - 1] > window_TH) {
                                    //pos.erase(pos.begin()-1);
                                    e_pos.push_back(g_pos[it - 1]);
                                }
                            } else {
                                e_pos.push_back(g_pos[it - 1]);
                            }
                            if (!(e_pos.empty())) {
                                if (g_pos[it + 1] - e_pos[e_pos.size() - 1] > window_TH) {
                                    //pos.erase(pos.begin()-1);
                                    e_pos.push_back(g_pos[it + 1]);
                                }
                            } else {
                                e_pos.push_back(g_pos[it + 1]);
                            }

                        } else {
                            if (!(e_pos.empty())) {
                                if (g_pos[it + 1] - e_pos[e_pos.size() - 1] > window_TH) {
                                    //pos.erase(pos.begin()-1);
                                    e_pos.push_back(g_pos[it + 1]);
                                }
                            } else {
                                e_pos.push_back(g_pos[it + 1]);
                            }

                        }
                        //it-1 and it+1 
                    }

                    it++;


                }
            } else {

                if (!(e_pos.empty())) {
                    if (g_pos[it - 1] - e_pos[e_pos.size() - 1] > window_TH) {
                        //pos.erase(pos.begin()-1);
                        e_pos.push_back(g_pos[it - 1]);
                    }
                } else {
                    e_pos.push_back(g_pos[it - 1]);
                }

                if (!(e_pos.empty())) {
                    if (g_pos[it] - e_pos[e_pos.size() - 1] > window_TH) {
                        //pos.erase(pos.begin()-1);
                        e_pos.push_back(g_pos[it]);
                    }
                } else {
                    e_pos.push_back(g_pos[it]);
                }


            }
            size++; // go to check last option when enter in loop (size>3) or as previous
        }
        if (size > 0 || size1 > 0) {
            if (g_pos[g_pos.size() - 1] - e_pos[e_pos.size() - 1] > window_TH) {
                //pos.erase(pos.begin()-1);
                e_pos.push_back(g_pos[g_pos.size() - 1]);
            } else {
                if (coverage[g_pos[g_pos.size() - 1] - 1] < coverage[e_pos[e_pos.size() - 1] - 1])
                    e_pos[e_pos.size() - 1] = g_pos[g_pos.size() - 1];
            }
        }
    }

    for (int it = 0; it < e_pos.size(); it++) {
        if (e_pos[it] >= st_end_th && e_pos[it] <= coverage.size() - st_end_th && coverage[e_pos[it]] <= cov_th) {
            f_pos.push_back(e_pos[it]);
        }
    }
    // cout << e_pos.size() << endl;

    return f_pos;

}

//11/11/19: not required for now

vector<int> mis_assembly_chimera::compare_singleSC(vector<int> singleSC_list, list<int> ex_bo, list<int> sc_pat) {
    // take single sc list and filter on the basis of exon-intron junction boundaries
    utils utils;
    vector<int> single_SC;
    single_SC.reserve(max_coverage);
    int boundary_margin = 5;
    
    int dip = 0, intron_counter = 1;
    for (int sc = 0; sc < singleSC_list.size(); sc++) {//no SCs in intronic region, so no check
        //  outfile << contig_name << "\t" << single_SC[sc] << "\t" << "By softclips" << endl;
        for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) {
            if ((singleSC_list[sc] >= (*boundary - boundary_margin)) && (singleSC_list[sc] <= (*boundary + boundary_margin))) {
                singleSC_list.erase(singleSC_list.begin() + sc); //pos< intron_start-100 or pos>intron_end+100
                sc = sc - 1;
                dip++;
                goto out_of_boundary_check;
            }

        }
        // *update 11/11/19*: not required because sc_pat contains all type of SC patterns now including single sc
        //        for (list<int>::iterator sc_it = sc_pat.begin(); sc_it != sc_pat.end(); sc_it++) { //relaxation window is 50 because for gradual check window size is 100, so to dilute the effect coverage change by exon boundary on window                  
        //            if ((singleSC_list[sc] >= (*sc_it - 5)) && (singleSC_list[sc] <= (*sc_it + 5))) {
        //                singleSC_list.erase(singleSC_list.begin() + sc); //pos< intron_start-100 or pos>intron_end+100
        //                sc = sc - 1;
        //                dip++;
        //                goto out_of_boundary_check;
        //            }
        //        }
out_of_boundary_check:
        ;
    }
    if (!singleSC_list.empty()) {
        for (int sc = 0; sc < singleSC_list.size(); sc++) {//no SCs in intronic region, so no check
            single_SC.push_back(singleSC_list[sc]);
        }
    }


    return single_SC;
}

vector<int> mis_assembly_chimera::compare_gradual_cov_change(vector<int>gradual_dip_chimera, list<int> ex_bo, vector<int> sc_pat) {
    utils utils;
    int dip = 0, intron_counter = 1;
    int boundary_margin = 5;
    //vector<int> gradual_cov_change;
  //  gradual_cov_change.reserve();
    //                                 
    //                                        }
    for (int it_window = 0; it_window < gradual_dip_chimera.size(); it_window++) {
        dip = 0;
        // to filter for splice-junction
        for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { //relaxation window is 50 because for gradual check window size is 100, so to dilute the effect coverage change by exon boundary on window                  
            if ((gradual_dip_chimera[it_window] >= (*boundary - boundary_margin)) && (gradual_dip_chimera[it_window] <= (*boundary + boundary_margin))) {
                gradual_dip_chimera.erase(gradual_dip_chimera.begin() + it_window);
                it_window = it_window - 1;
                intron_counter = 1;
                dip++;
                goto outer_dip;
            }
        }
        // to filter for miss-assembled regions
        for (vector<int>::iterator sc = sc_pat.begin(); sc != sc_pat.end(); sc++) { //relaxation window is 50 because for gradual check window size is 100, so to dilute the effect coverage change by exon boundary on window                  
            if (abs(gradual_dip_chimera[it_window] - *sc) <= consecutive_missAssembled_pos_dist) {
                gradual_dip_chimera.erase(gradual_dip_chimera.begin() + it_window);
                it_window = it_window - 1;
                intron_counter = 1;
                dip++;
                goto outer_dip;
            }
        }
outer_dip:
        ;
    }
    //return gradual_cov_change; //its not even getting any value? 16 dec2020
    return gradual_dip_chimera;

}
// list<int> mis_assembly_chimera::crossed_sc_pattern(string contig_name, vector<int>left_SC_base_pos, vector<string>left_SC_seq, vector<int>right_SC_base_pos, vector<string>right_SC_seq, vector<int>cov) {
///*extract crossing soft-clipped patterns b/c considering distance and overlap of 10 bases.
// Because positions for sofctlips pattersn facing each other has been cater by hisat alignment in the form of slice-junctions*/
//    
//     cout<<"entering sc pattern"<<endl;
//    std::map<int, int>::iterator itr1_cov;
//    int last_leftSC = 0, last_rightSC = 0;
//    float consective_leftSC = 0, consective_rightSC = 0;
//    float previous_last_leftSC = 0, previous_last_rightSC = 0;
//    vector<string> left_seq, right_seq;
//    list<int> SC_chimera;
//    // cout<< contig_name<<" in soft clip1 "<<cov.size()<<" "<<endl;
//    // cout<<"entering SC"<<endl;
//
//    // try{
//    for (int m = 1; m < left_SC_base_pos.size(); m++) {
//        if (left_SC_base_pos[m] - left_SC_base_pos[m - 1] >= 5) //if at a distance
//        {
//            //   cout<<"there is a diff"<<endl;
//            left_SC_base_pos[m - 1] = 0;
//            left_SC_seq[m - 1] = "";
//        } else {
//            continue;
//        }
//    }
//    for (int k = 1; k < right_SC_base_pos.size(); k++) {
//        if (right_SC_base_pos[k] - right_SC_base_pos[k - 1 ] >= 5) { //5?
//            right_SC_base_pos[k - 1] = 0;
//            right_SC_seq[k - 1] = " ";
//        } else {
//            continue;
//        }
//    }
//    for (int mm = 1; mm < left_SC_base_pos.size(); mm++) { //for 1 left SC cluster check all right soft clip clusters ...between two clusters 0 is present to differentiate them
//        if (left_SC_base_pos[mm] != 0) {
//            last_leftSC = left_SC_base_pos[mm];
//            consective_leftSC += 1;
//            left_seq.push_back(left_SC_seq[mm]);
//        } else {
//            if (consective_leftSC > 4) {// && max_consec_leftSC < consective_leftSC) {
//
//                for (int kk = 1; kk < right_SC_base_pos.size(); kk++) {
//                    //   cout<<right_SC_base_pos[kk]<<endl;
//                    if (right_SC_base_pos[kk] != 0) {
//                        last_rightSC = right_SC_base_pos[kk];
//                        consective_rightSC += 1;
//                        right_seq.push_back(right_SC_seq[kk]);
//
//                    } else {
//                        if (consective_rightSC > 4)// && max_consec_rightSC < consective_rightSC) //max cluster supported by  >4 SC reads
//                        {
//                            if (previous_last_leftSC != last_leftSC && previous_last_rightSC != last_rightSC) { //don't repeat same base_positin
//                                if (last_leftSC - last_rightSC <= 5 && last_leftSC - last_rightSC >= -5) {//DISTANCE AND OVErlap    
//                                    previous_last_rightSC = last_rightSC;
//                                    previous_last_leftSC = last_leftSC;
//                                    for (vector<int>::iterator c = cov.begin() ; c != cov.end(); c++)  {
//                                        int i = std::distance(cov.begin(), c) - 1; // 1 because cov index from 0;
//                                        if ((i == previous_last_leftSC && (consective_leftSC / *c)*100 >= SC_support_TH) || (i == previous_last_rightSC && (consective_rightSC / (*c) * 100 >= SC_support_TH))) {
//                                            SC_chimera.push_back(previous_last_rightSC);
//                                            SC_chimera.push_back(previous_last_leftSC);
//                                        }
//
//                                    }
//                                }
//                            }
//                        }
//                        right_seq.clear();
//                        last_rightSC = 0;
//                        consective_rightSC = 0;
//                    }
//                }
//
//            }
//            left_seq.clear();
//            last_leftSC = 0;
//            consective_leftSC = 0;
//        }
//    }
//    // cout<<"leaving SC"<<endl;
//    return SC_chimera;
//}

vector<int> mis_assembly_chimera::filter_singleSC(vector<int> single_SC, map <int, string> unmappedMate_reads, map<int, string> distMapped_reads) {
    /*filter single sided softclips with unmapped reads and distantly mapped reads arounds sc positions, these patterns will appear only in hisat alignment*/
    utils utils;
    //cout << "entering sc pattern" << endl;
    std::map<int, int>::iterator itr1_cov;
    int last_leftSC = 0, last_unmapped = 0;
    float consective_leftSC = 0, unmapped_count = 1, dist_mapped_count = 1;
    float previous_last_leftSC = 0, previous_last_rightSC = 0;
    
 
    vector<string> unmapped_dir, dist_mapped_dir;
    unmapped_dir.reserve(50); dist_mapped_dir.reserve(50);
    
    vector<int> SC_chimera, unmapped_pos, dist_mapped_pos;
    SC_chimera.reserve(500); unmapped_pos.reserve(50); dist_mapped_pos.reserve(50);
    
    // cout<< contig_name<<" in soft clip1 "<<cov.size()<<" "<<endl;
    // cout<<"entering SC"<<endl;

    // try{
    for (map <int, string>::iterator it = unmappedMate_reads.begin(); it != unmappedMate_reads.end(); it++) {
        unmapped_pos.push_back(it->first);
        unmapped_dir.push_back(it->second);
    }

    for (map <int, string>::iterator itt = distMapped_reads.begin(); itt != distMapped_reads.end(); itt++) {
        dist_mapped_pos.push_back(itt->first);
        dist_mapped_dir.push_back(itt->second);
    }

    for (int m = 1; m < unmapped_pos.size(); m++) {
        if (unmapped_pos[m] - unmapped_pos[m - 1] >= 25) //if at a distance
        {
            unmapped_pos[m - 1] = 0;
            unmapped_dir[m - 1] = " ";
        } else {
            continue;
        }
    }
    for (int k = 1; k < dist_mapped_pos.size(); k++) {
        if (dist_mapped_pos[k] - dist_mapped_pos[k - 1 ] >= 25) { //5?
            dist_mapped_pos[k - 1] = 0;
            dist_mapped_dir[k - 1] = " ";
        } else {
            continue;
        }
    }
    for (int n = 0; n < single_SC.size(); n++) { //for 1 left SC cluster check all right soft clip clusters ...between two clusters 0 is present to differentiate them
        for (int kk = 0; kk < unmapped_pos.size(); kk++) {
            //   cout<<right_SC_base_pos[kk]<<endl;
            if (unmapped_pos[kk] != 0) {
                if (unmapped_pos[kk] - single_SC[n] <= 50 && unmapped_pos[kk] - single_SC[n] >= 0 && unmapped_dir [kk] == "R") {//DISTANCE AND OVErlap    // half of the read length 
                    unmapped_count++;
                } else if (unmapped_pos[kk] - single_SC[n] > -150 && unmapped_pos[kk] - single_SC[n] <= 0 && unmapped_dir [kk] == "F") {//DISTANCE AND OVErlap    
                    unmapped_count++;
                }
            }
        }

        for (int mm = 0; mm < dist_mapped_pos.size(); mm++) {
            //   cout<<right_SC_base_pos[kk]<<endl;
            if (dist_mapped_pos[mm] != 0) {
                if (dist_mapped_pos[mm] - single_SC[n] <= 50 && dist_mapped_pos[mm] - single_SC[n] >= 0 && dist_mapped_dir [mm] == "R") {//DISTANCE AND OVErlap    
                    dist_mapped_count++;
                } else if (dist_mapped_pos[mm] - single_SC[n] > -150 && dist_mapped_pos[mm] - single_SC[n] <= 0 && dist_mapped_dir [mm] == "F") {//DISTANCE AND OVErlap // read start pos + read length ~100 + distance    
                    dist_mapped_count++;
                }
            }
        }

        if (unmapped_count >= 3 || dist_mapped_count >= 3) //single sc pos surrounded by either unmapped reds or distantly mapped reads  to be sure about chimera.
            SC_chimera.push_back(single_SC[n]);
        unmapped_count = 0;
        dist_mapped_count = 0;
    }
    // cout<<"leaving SC"<<endl;
    return SC_chimera;
}

vector<int> mis_assembly_chimera::sc_pattern(string contig_name, string contig_seq, multimap<int, string> left_SC_pos_seq, multimap<int, string> right_SC_pos_seq, vector<int>cov, string &output_file, string PATH, string bam_file) {
    /*1. extract crossed soft-clipped patterns b/c considering distance and overlap of 5 bases.
      2. extract SCs pattern towards each other if their seq matched on opposite sites in contig
      3. extract single-sided softclips 
     * return map with pattern type as key and vector for respective positions */

    float SC_support_TH2 = float (SC_support_TH / 2);// +  35; // to make it 40 for this module
    utils utils;
    
    //int SC_support_TH = 75;
   //float SC_support_TH2 = SC_support_TH -  35;
  // float SC_support_TH3 = SC_support_TH -  15; //60
    std::replace(cov.begin(), cov.end(), 0, 1); // replace zero cov to 1 to avoid floating point exception
    //cout << "entering sc pattern" << endl;
    std::map<int, int>::iterator itr1_cov;
    vector <int>::iterator it;
    int last_leftSC = 0, last_rightSC = 0, found_right, found_right_end, found_left, found_left_end, a = 0, left_st_ind, left_en_ind;
    float consective_leftSC = 0, consective_rightSC = 0;
    int consective_SC_TH = 3;
    
    float previous_last_leftSC = 0, previous_last_rightSC = 0;
    string right_cons_seq, left_cons_seq, blast_found;
    int cons_seq_size_TH = 10;
    
    string first_frag, sec_frag, query_frag, strand;
    int query_frag_avg_TH = 25;
    int left_right_dist = 10;
    int repeat_TH = 100;
    
    vector<string> left_seq, right_seq, left_SC_seq, right_SC_seq;
    left_seq.reserve(max_coverage), right_seq.reserve(max_coverage), left_SC_seq.reserve(max_coverage), right_SC_seq.reserve(max_coverage);
    
    vector<int> sc_patterns, left_SC_base_pos, right_SC_base_pos;
    sc_patterns.reserve(max_coverage), left_SC_base_pos.reserve(max_coverage); right_SC_base_pos.reserve(max_coverage);
   
    bool isMatchFound, problem_identified = false, towards_found = false, crossed_found = false, single_sc = false, toleft = true;
    
    ofstream out_file;
    out_file.open(output_file.c_str(), fstream::app);
    // cout<<"entering SC"<<endl;
    // try{
    //cout<< contig_name<<endl;
    for (multimap<int, string>::iterator Values = left_SC_pos_seq.begin(); Values != left_SC_pos_seq.end(); ++Values) {
        left_SC_base_pos.push_back(Values->first);
        left_SC_seq.push_back(Values->second);
        //            cout << (*Values).first << " is married to " <<
        //                (*Values).second << endl;
    }

    for (multimap<int, string>::iterator Values = right_SC_pos_seq.begin(); Values != right_SC_pos_seq.end(); ++Values) {
        right_SC_base_pos.push_back(Values->first);
        right_SC_seq.push_back(Values->second);
    }

    for (int m = 1; m < left_SC_base_pos.size(); m++) {
        if (left_SC_base_pos[m] != left_SC_base_pos[m - 1]) //if at a distance
        {
            left_SC_base_pos.insert(left_SC_base_pos.begin() + m, 0);
            left_SC_seq.insert(left_SC_seq.begin() + m, "");
            m++;
        } else {
            continue;
        }
    }
    left_SC_base_pos.push_back(0); // ti deal with last SC pattern
    for (int k = 1; k < right_SC_base_pos.size(); k++) {
        if (right_SC_base_pos[k] != right_SC_base_pos[k - 1 ]) { //5?
            right_SC_base_pos.insert(right_SC_base_pos.begin() + k, 0);
            right_SC_seq.insert(right_SC_seq.begin() + k, "");
            k++;
        } else {
            continue;
        }
    }
    right_SC_base_pos.push_back(0); // it deal with last SC pattern
    
    if (!(right_SC_base_pos.empty())) {
        
        for (int mm = 0; mm < right_SC_base_pos.size(); mm++) { //for 1 right SC cluster check all left soft clip clusters ...between two clusters 0 is present to differentiate them
           
            if (right_SC_base_pos[mm] != 0) {
                last_rightSC = right_SC_base_pos[mm];
                consective_rightSC += 1;
                right_seq.push_back(right_SC_seq[mm]);
          
            } else {
                if (consective_rightSC >= consective_SC_TH && abs(last_rightSC - previous_last_rightSC) > consecutive_missAssembled_pos_dist) {
                    
                    int pos_check;
                   
                    if (cov.at(last_rightSC - 1) > 0) // to avoid zero cov at sc position
                        pos_check = last_rightSC ;
                    else 
                        pos_check = last_rightSC - 1; // -2 because its RSC
                     //cout << pos_check << " " << cov.at(pos_check - 1) << " " << consective_rightSC / (cov.at(pos_check - 1)) * 100 << endl;
                    
                    if (float(consective_rightSC / cov.at(pos_check - 1)) * 100 >= SC_support_TH2) { // -1 because coverage starts from 0 index not one
                        right_cons_seq = utils.consensus_seq_CAP3(right_seq, !toleft, PATH);
                        
                        if (right_cons_seq.size() >= cons_seq_size_TH) { // because less than done blast doesn't work
                            previous_last_rightSC = last_rightSC;
                            //
                            string blast_hit = perform_blast(right_cons_seq, contig_seq); //boost::regex_search(contig_seq, present, regexPattern);
                            vector<string> temp;
                            utils.str_split(blast_hit, temp, delimiter);
                            found_right = atoi(temp[0].c_str());
                            //found_right_end = atoi(temp[2].c_str()) - found_right;
                            strand = temp[1];

                            if (found_right != -1) {
                                if (strand == "+") {
                                    //found_right = present.position(a); // for right start pos of match should be same as start of LSC (start of read containing LSC)
                                    found_right_end = atoi(temp[2].c_str()) - found_right;
                                    // check for towards pattern if RSc seq found somewhere

                                    for (int kk = 0; kk < left_SC_base_pos.size(); kk++) {
                                        //   cout<<right_SC_base_pos[kk]<<endl;
                                        if (left_SC_base_pos[kk] != 0) {
                                            last_leftSC = left_SC_base_pos[kk];
                                            consective_leftSC += 1;
                                            left_seq.push_back(left_SC_seq[kk]);

                                        } else {
                                           
                                            if (consective_leftSC >= consective_SC_TH && previous_last_leftSC != last_leftSC) {// && max_consec_rightSC < consective_rightSC) //max cluster supported by  >4 SC reads
                                                int pos_check;
                                                
                                                if (cov.at(last_leftSC - 1 ) > 1) // to avoid zero cov at sc position
                                                    pos_check = last_leftSC ;
                                                else 
                                                    pos_check = last_leftSC - 1;
                                                
                                                if (float(consective_leftSC / cov.at(pos_check - 1)) * 100 >= SC_support_TH2) { // -1 because coverage starts from 0 index not one
                                                    left_st_ind = kk - consective_leftSC;
                                                    left_en_ind = kk - 1;
                                                    left_cons_seq = utils.consensus_seq_CAP3(left_seq, toleft, PATH); // "left");
                                                    
                                                    if (left_cons_seq.size() >= cons_seq_size_TH) {
                                                        previous_last_leftSC = last_leftSC;

                                                        string blast_hit_left = perform_blast(left_cons_seq, contig_seq); // boost::regex_search(contig_seq, present, regexPattern);
                                                        vector<string> temp2;
                                                        utils.str_split(blast_hit_left, temp2, delimiter);
                                                        found_left = atoi(temp2[0].c_str());
                                                        found_left_end = atoi(temp2[2].c_str());// 17jan - found_left;
                                                        string strand_left = temp2[1];
                                                        
                                                        if (found_left != -1 && abs(found_left_end - last_leftSC) > left_right_dist && strand_left == "+" ){ //if (found_left != -1 && strand_left == "+") { // for +ve strand of RSC only +ve LSC will be searched
                                                            // if(strand == "+")
                                                            found_left = found_left; // for left end pos of match should be same as RSC start pos
                                                           
                                                            if ((abs(found_right - last_leftSC) <= left_right_dist) || (abs(found_left - last_rightSC) <= left_right_dist)) { //towards pattern, seq match may not at exact LSC start 
                                                                
                                                                int start_region = last_rightSC;
                                                                int end_region = found_right;
                                                                //*******************************************************************************************
                                                                query_frag = contig_seq.substr(start_region, end_region - start_region); // middle suspicious fragment
                                                                first_frag = contig_seq.substr(0, start_region - 1); // fragment before suspicious fragment
                                                                
                                                                if (end_region + 1 < contig_seq.size())
                                                                    sec_frag = contig_seq.substr(end_region + 1, contig_seq.size() - 1); // fragment after
                                                                else
                                                                    sec_frag = "";
                                                                //sec_frag = contig_seq.substr(end_region + 1, contig_seq.size()); // fragment after
                                                                string target = first_frag + string(query_frag.size(), 'N') + sec_frag;
                                                                string blast_hit = perform_blast(query_frag, target, " ");
                                                                
                                                                vector<string> temp2; // = utils.split(blast_hit, "\t");
                                                                utils.str_split(blast_hit, temp2, delimiter);
                                                                blast_found = temp2[0];
                                                                int target_start = atoi(temp2[1].c_str());
                                                                int target_end = atoi(temp2[2].c_str());
                                                                bool repeat = false;
                                                                
                                                                if (!sc_patterns.empty()) { // check repeat positions
                                                                   
                                                                    for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                      
                                                                        if (abs(target_start - *sc) < repeat_TH) {
                                                                          
                                                                            for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                              
                                                                                if ((abs(target_end - *sc) < repeat_TH)) {
                                                                                    repeat = true;
                                                                                    break;
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                                if (repeat) { // skip if fragment is already dealt (to avoid repetition of self chimeric fragments)
                                                                    //                                                                    sc_patterns.push_back(target_start);
                                                                    //                                                                    sc_patterns.push_back(target_end);
                                                                    goto out_right_next;
                                                                } else if (blast_found != "*") { //  hit found
                                                                    sc_patterns.push_back(target_start);
                                                                    sc_patterns.push_back(target_end);
                                                                    
                                                                    if (blast_found == "+") {
                                                                        problem_identified = true; //selfChimera
                                                                        sc_patterns.push_back(start_region);
                                                                        sc_patterns.push_back(end_region);
                                                                        
                                                                        if ((abs(target_start - start_region) > cons_seq_size_TH) && abs(target_end - end_region) > cons_seq_size_TH) {
                                                                            
                                                                            out_file << contig_name << "\t" << "selfChimera" << "\t" << start_region << "-" << end_region << "\t" << "forward" << "\t" << "discard" << endl;
                                                                            // erase treated left sc, so that won't check them later
                                                                            left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                            left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                            kk = kk - consective_leftSC;
                                                                            right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                            right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                            mm = mm - consective_rightSC;
                                                                        }
                                                                        goto out_right_next;
                                                                    } else { // RC of query seq strand -
                                                                        query_frag = utils.Rcomplement(query_frag);
                                                                        //blast_found = perform_blast(query_frag, first_frag, sec_frag);
                                                                        problem_identified = true; //selfChimera   
                                                                        sc_patterns.push_back(start_region);
                                                                        sc_patterns.push_back(end_region);
                                                                        out_file << contig_name << "\t" << "selfChimera" << "\t" << start_region << "-" << end_region << "\t" << "Reverse complement" << "\t" << "discard" << endl;
                                                                        //id type pos (st-end) RC
                                                                        //erase treated left sc, so that won't check them later
                                                                        left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                        left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                        kk = kk - consective_leftSC; // update index for LSC
                                                                        right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                        right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                        mm = mm - consective_rightSC;
                                                                        goto out_right_next;
                                                                    }
                                                                } else { // NO blast hit found so check any SC pattern initiate in emitted region
                                                                    // left_found = lsc, right_found = rsc
                                                                    float lsc_occurance = 0;
                                                                    //it = find(left_SC_base_pos.begin(), left_SC_base_pos.end(), start_region);
                                                                    bool left_found_in_right = false;

                                                                    for (it = left_SC_base_pos.begin(); it != left_SC_base_pos.end(); it++) {
                                                                        if (abs(*it - start_region) <= left_right_dist) {
                                                                            left_found_in_right = true;
                                                                            break;
                                                                        }
                                                                    }

                                                                    if(left_found_in_right) { //if (it != left_SC_base_pos.end()) {
                                                                        int pos = std::distance(left_SC_base_pos.begin(), it);
                                                                        vector <string> seq_lsc_occ;
                                                                        
                                                                        for (int l_found = pos; left_SC_base_pos[l_found] != 0; l_found++) {
                                                                            seq_lsc_occ.push_back(left_SC_seq[l_found]);
                                                                            lsc_occurance++;
                                                                        }
                                                                        if (lsc_occurance >= consective_SC_TH && previous_last_leftSC != end_region) { // lsc pattern found at the start of emitted sequence
                                                                            int pos_check;
                                                                            
                                                                            if (cov.at(start_region ) > 1) // to avoid zero cov at sc position
                                                                                pos_check = start_region ;
                                                                            else 
                                                                                pos_check = start_region - 1;
                                                                            
                                                                          //  cout << cov.at(pos_check - 1) <<  "  " << lsc_occurance << " "<< float(lsc_occurance / cov.at(pos_check - 1)) * 100 << endl;
                                                                            
                                                                            if (float(lsc_occurance / cov.at(pos_check - 1)) * 100 >= SC_support_TH2) { // -1 because coverage starts from 0 index not one // check cov for SC found at start_region
                                                                              
                                                                                string cons_lsc_occ = utils.consensus_seq_CAP3(seq_lsc_occ, toleft, PATH); //"left");
                                                                                
                                                                                if (cons_lsc_occ.size() >= cons_seq_size_TH) {

                                                                                    string left_found_blast = perform_blast(cons_lsc_occ, contig_seq); //boost::regex_search(contig_seq, present, regexPattern);
                                                                                    vector<string> temp3; // = utils.split(left_found_blast, "\t");
                                                                                    utils.str_split(left_found_blast, temp3, delimiter);
                                                                                    int cons_left_found = atoi(temp3[0].c_str());
                                                                                    int cons_left_end = atoi(temp3[2].c_str()) - cons_left_found;
                                                                                    string strand_left_found = temp3[1];
                                                                                    
                                                                                    if (cons_left_found != -1) {
                                                                                        
                                                                                        if (strand_left_found == "+") {
                                                                                            
                                                                                            cons_left_found = cons_left_found + cons_left_end; // end pos of match to insert
                                                                                            problem_identified = true; // *translocation*
                                                                                            sc_patterns.push_back(start_region);
                                                                                            sc_patterns.push_back(end_region);
                                                                                            sc_patterns.push_back(cons_left_found);
                                                                                            out_file << contig_name << "\t" << "translocation" << "\t" << start_region << "-" << end_region << "\t" << "forward" << "\t" << "insert at:" << cons_left_found << endl;
                                                                                            // erase treated left sc, so that won't check them later
                                                                                            left_st_ind = pos;
                                                                                            left_en_ind = pos + lsc_occurance;
                                                                                            left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                                            left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                                            right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                                            right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                                            mm = mm - consective_rightSC;
                                                                                            goto out_right_next;
                                                                                            
                                                                                        } else{// if (strand_left_found == "-") {
                                                                                            cons_lsc_occ = utils.Rcomplement(cons_lsc_occ); // take reverse complement of LSC seq 

                                                                                            //int cons_left_found = perform_blast(cons_lsc_occ, contig_seq); //boost::regex_search(contig_seq, present, regexPattern);
                                                                                            //int cons_left_found = present.position(a); // start pos of match to insert because RC
                                                                                            problem_identified = true; // *inversion and translocation*
                                                                                            sc_patterns.push_back(start_region);
                                                                                            sc_patterns.push_back(end_region);
                                                                                            sc_patterns.push_back(cons_left_found);
                                                                                            out_file << contig_name << "\t" << "inversion_translocation" << "\t" << start_region << "-" << end_region << "\t" << "Reverse complement" << "\t" << "insert at:" << cons_left_found << endl;
                                                                                            // erase treated lesft sc, so that won't check them later
                                                                                            left_st_ind = pos;
                                                                                            left_en_ind = pos + lsc_occurance;
                                                                                            left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                                            left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                                            right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                                            right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                                            mm = mm - consective_rightSC;
                                                                                            goto out_right_next;

                                                                                        }
                                                                                    }
                                                                                } // consense check
                                                                            } // cov check
                                                                        }
                                                                        lsc_occurance = 0;
                                                                    }
                                                                    if (!problem_identified) { // no LSC at start found so now check right corner of emitted sequence have any RSC pattern
                                                                        float rsc_occurance = 0;
                                                                        bool right_found_in_right = false;
                                                                        
                                                                        for(it = right_SC_base_pos.begin(); it != right_SC_base_pos.end(); it++){
                                                                            
                                                                            if(abs(*it - found_right) <= left_right_dist){
                                                                                right_found_in_right = true;
                                                                                break;
                                                                            }
                                                                        }
                                                                        //it = find(right_SC_base_pos.begin(), right_SC_base_pos.end(), found_right);
                                                                        if (right_found_in_right){ //if (it != right_SC_base_pos.end()) {
                                                                            
                                                                            int pos = std::distance(right_SC_base_pos.begin(), it);
                                                                            vector <string> seq_rsc_occ;
                                                                            
                                                                            for (int r_found = pos; right_SC_base_pos[r_found] != 0; r_found++) {
                                                                                seq_rsc_occ.push_back(right_SC_seq[r_found]);
                                                                                rsc_occurance++;
                                                                            }
                                                                            if (rsc_occurance >= consective_SC_TH && found_right != previous_last_rightSC) { // rsc pattern found at the end of emitted sequence
                                                                                int pos_check;
                                                                               
                                                                                if (cov.at(found_right - 1) > 1) // to avoid zero cov at sc position
                                                                                    pos_check = found_right ;
                                                                                else 
                                                                                    pos_check = found_right - 1;

                                                                                if (float(rsc_occurance / cov.at(pos_check - 1)) * 100 >= SC_support_TH2) { // -1 because coverage starts from 0 index not one // outward sc found at position where RSC seq found in contig 
                                                                                    
                                                                                    string cons_rsc_occ = utils.consensus_seq_CAP3(seq_rsc_occ, !toleft, PATH) ; //"right");

                                                                                    if (cons_rsc_occ.size() >= cons_seq_size_TH) {
                                                                                        
                                                                                        string cons_rsc_blast = perform_blast(cons_rsc_occ, contig_seq); //boost::regex_search(contig_seq, present, regexPattern);
                                                                                        vector<string> temp3; // = utils.split(cons_rsc_blast, "\t");
                                                                                        utils.str_split(cons_rsc_blast, temp3, delimiter);
                                                                                        int cons_right_found = atoi(temp3[0].c_str());
                                                                                        int cons_right_end = atoi(temp3[2].c_str()) - cons_right_found;
                                                                                        string cons_right_strand = temp3[1];

                                                                                        if (cons_right_found != -1) {
                                                                                            
                                                                                            if (cons_right_strand == "+") {
                                                                                                //int cons_right_found = present.position(a); //start pos of match to insert because RC
                                                                                                problem_identified = true; // *translocation*
                                                                                                sc_patterns.push_back(start_region);
                                                                                                sc_patterns.push_back(end_region);
                                                                                                sc_patterns.push_back(cons_right_found);
                                                                                                out_file << contig_name << "\t" << "translocation" << "\t" << start_region << "-" << end_region << "\t" << "forward" << "\t" << "insert at:" << cons_right_found << endl;
                                                                                                right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                                                right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                                                mm = mm - consective_rightSC;
                                                                                                goto out_right_next;
                                                                                                
                                                                                            } else{// if (cons_right_strand == "-") {
                                                                                                cons_rsc_occ = utils.Rcomplement(cons_rsc_occ); // take reverse complement of LSC seq 
                                                                                                cons_right_found = cons_right_found + cons_right_end ; // start pos of match to insert because RC
                                                                                                problem_identified = true; // *inversion and translocation*
                                                                                                sc_patterns.push_back(start_region);
                                                                                                sc_patterns.push_back(end_region);
                                                                                                sc_patterns.push_back(cons_right_found);
                                                                                                out_file << contig_name << "\t" << "inversion_translocation" << "\t" << start_region << "-" << end_region << "\t" << "Reverse complement" << "\t" << "insert at:" << cons_right_found << endl;
                                                                                                right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                                                right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                                                mm = mm - consective_rightSC;
                                                                                                goto out_right_next;
                                                                                            }
                                                                                        } //cov check
                                                                                    } //consenses check
                                                                                }
                                                                                rsc_occurance = 0;
                                                                            }
                                                                        }
                                                                    }
                                                                    //}
                                                                    //******************************************************************************
                                                                    if (!problem_identified) { // both RSC and LSC found but no sign of selfChimera, so now check for multigeneChimera before going towards single-sided SC, because it means only RSC found no LSC found(which we found in this case)

                                                                        int sum = 0;
//                                                                        for (int aa = start_region - 110; aa <= start_region - 10; aa++) {
//                                                                            sum = sum +  cov[aa];
//                                                                        }
//                                                                        int left_avg_cov = sum / 100;
//                                                                        sum = 0;
//
//                                                                        for (int aa = end_region + 10; aa <= end_region + 110; aa++) {
//                                                                            sum = sum + cov[aa];
//                                                                        }
//                                                                        int right_avg_cov = sum / 100;
//                                                                        sum = 0;
//
                                                                        for (int aa = start_region; aa <= end_region; aa++) {
                                                                            sum = sum + cov[aa];
                                                                        }
                                                                        int query_frag_avg = float(sum) / float(query_frag.size());

                                                                        if (abs(end_region - start_region) >= ignore_short_seq && query_frag_avg >= query_frag_avg_TH) { // read coverage should 25%  of length
                                                                            bool repeat = false;
                                                                            
                                                                            if (!sc_patterns.empty()) { // check miss-assembled regions
                                                                                
                                                                                for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                                    
                                                                                    if (abs(start_region - *sc) < repeat_TH) {
                                                                                        
                                                                                        for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                                            
                                                                                            if ((abs(end_region - *sc) < repeat_TH)) {
                                                                                                repeat = true;
                                                                                                break;
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                            if (repeat) {
                                                                                goto out_right_next;
                                                                            }

                                                                            //28oct21 avoid chimeras at the start or end positions of overlapped merged positions found in RI and merged by SCs
                                                                            float SC_TH_mc;
                                                                 
                                                                            if(find (left_SC_base_pos.begin(), left_SC_base_pos.end(), start_region - 1) != left_SC_base_pos.end()){ // if crisscross pattern
                                                                                SC_TH_mc = SC_support_TH2;
                                                                            }
                                                                            else // single sided pattern
                                                                                SC_TH_mc = SC_support_TH;
                                                                            //cout <<float(consective_rightSC / cov.at(pos_check - 1)) << endl;
                                                                            if (float(consective_rightSC / cov.at(pos_check - 1)) * 100 >= SC_TH_mc) { // 10feb22

                                                                            if (merged_contigs_itr == merged_contigs.end()) {
                                                                                problem_identified = true; //trasns chimera
                                                                                sc_patterns.push_back(start_region);
                                                                                sc_patterns.push_back(end_region);
                                                                                out_file << contig_name << "\t" << "multigeneChimera" << "\t" << start_region << "-" << end_region << "\t" << "forward" << "\t" << "new contig" << endl;
                                                                                // erase treated lesft sc, so that won't check them later
                                                                                left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                                left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                                kk = kk - consective_leftSC;
                                                                                right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                                right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                                mm = mm - consective_rightSC;

                                                                            } else { //identified chimera positions must be at one side of merged points in the merged contigs
                                                                                
                                                                                if ((start_region < start_overlap_merged && end_region < start_overlap_merged) || (start_region > end_overlap_merged && end_region > end_overlap_merged)) {
                                                                                    
                                                                                    problem_identified = true; //trasns chimera
                                                                                    sc_patterns.push_back(start_region);
                                                                                    sc_patterns.push_back(end_region);
                                                                                    out_file << contig_name << "\t" << "multigeneChimera" << "\t" << start_region << "-" << end_region << "\t" << "forward" << "\t" << "new contig" << endl;
                                                                                    // erase treated lesft sc, so that won't check them later
                                                                                    left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                                    left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                                    kk = kk - consective_leftSC;
                                                                                    right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                                    right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                                    mm = mm - consective_rightSC;
                                                                                } // else go out
                                                                            }
                                                                            }
                                                                            goto out_right_next;

                                                                        } else if (query_frag.size() < ignore_short_seq) {

                                                                                problem_identified = true;
                                                                                sc_patterns.push_back(start_region);
                                                                            sc_patterns.push_back(end_region);
                                                                            out_file << contig_name << "\t" << "unsupported_insertion" << "\t" << start_region << "-" << end_region << "\t" << "forward" << "\t" << "Discard" << endl;
                                                                            left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                            left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                            kk = kk - consective_leftSC;

                                                                            right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                            right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                            mm = mm - consective_rightSC;
                                                                            goto out_right_next;
                                                                        }

                                                                    }
                                                                }
                                                                //added 1/12
                                                                if (!problem_identified) { // both RSC and LSC found but no sign of selfChimera, so now check for multigeneChimera before going towards single-sided SC, because it means only RSC found no LSC found(which we found in this case)

                                                                    int sum =  0;
//                                                                    for (int aa = start_region - 110; aa <= start_region - 10; aa++) {
//                                                                        sum = sum + cov[aa];
//                                                                    }
//                                                                    float left_avg_cov = sum / 100;
//                                                                    sum = 0;
//
//                                                                    for (int aa = end_region + 10; aa <= end_region + 110; aa++) {
//                                                                        sum = sum + cov[aa];
//                                                                    }
//                                                                    float right_avg_cov = sum / 100;
//                                                                    sum = 0;

                                                                    for (int aa = start_region; aa <= end_region; aa++) {
                                                                        sum = sum + cov[aa];
                                                                    }
                                                                    float query_frag_avg = float(sum) / float(query_frag.size());

                                                                    if (abs(end_region - start_region) >= ignore_short_seq && query_frag_avg >= query_frag_avg_TH) {
                                                                        bool repeat = false;
                                                                        
                                                                        if (!sc_patterns.empty()) { // check miss-assembled regions
                                                                            
                                                                            for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                                
                                                                                if (abs(start_region - *sc) < repeat_TH) {
                                                                                    
                                                                                    for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                                        
                                                                                        if ((abs(end_region - *sc) < repeat_TH)) {
                                                                                            repeat = true;
                                                                                            break;
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                        if (repeat) {
                                                                            goto out_right_next;
                                                                        }
                                                                        //28oct21 avoid chimeras at the start or end positions of overlapped merged positions found in RI and merged by SCs                                                                        
                                                                        float SC_TH_mc;

                                                                        if (find(left_SC_base_pos.begin(), left_SC_base_pos.end(), start_region - 1) != left_SC_base_pos.end()) { // if crisscross pattern
                                                                            SC_TH_mc = SC_support_TH2;
                                                                        } else // single sided pattern
                                                                            SC_TH_mc = SC_support_TH;
                                                                        //cout <<float(consective_rightSC / cov.at(pos_check - 1)) << endl;           
                                                                        if (float(consective_rightSC / cov.at(pos_check - 1)) * 100 >= SC_TH_mc) { // 10feb22
                                                                            if (merged_contigs_itr == merged_contigs.end()) {

                                                                                problem_identified = true; //trasns chimera 
                                                                                sc_patterns.push_back(start_region);
                                                                                sc_patterns.push_back(end_region);
                                                                                out_file << contig_name << "\t" << "multigeneChimera" << "\t" << start_region << "-" << end_region << "\t" << "forward" << "\t" << "new contig" << endl;
                                                                                // erase treated lesft sc, so that won't check them later
                                                                                left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                                left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                                kk = kk - consective_leftSC;
                                                                            right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                            right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                            mm = mm - consective_rightSC;

                                                                        } else { //identified chimera positions must be at one side of merged points in the merged contigs

                                                                            if ((start_region < start_overlap_merged && end_region < start_overlap_merged) || (start_region > end_overlap_merged && end_region > end_overlap_merged)) {
                                                                                
                                                                                problem_identified = true; //trasns chimera 
                                                                                sc_patterns.push_back(start_region);
                                                                                sc_patterns.push_back(end_region);
                                                                                out_file << contig_name << "\t" << "multigeneChimera" << "\t" << start_region << "-" << end_region << "\t" << "forward" << "\t" << "new contig" << endl;
                                                                                // erase treated lesft sc, so that won't check them later
                                                                                left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                                left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                                kk = kk - consective_leftSC;
                                                                                right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                                right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                                mm = mm - consective_rightSC;
                                                                            } // else go out
                                                                        }
                                                                    }
                                                                        goto out_right_next;
                                                                        
                                                                    } else if (query_frag.size() < ignore_short_seq) {
                                                                        problem_identified = true;
                                                                        sc_patterns.push_back(start_region);
                                                                        sc_patterns.push_back(end_region);
                                                                        out_file << contig_name << "\t" << "unsupported_insertion" << "\t" << start_region << "-" << end_region << "\t" << "forward" << "\t" << "Discard" << endl;
                                                                        left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                        left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                        kk = kk - consective_leftSC;

                                                                        right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                        right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                        mm = mm - consective_rightSC;
                                                                        goto out_right_next;
                                                                    }
                                                                }

                                                            } // not facing each other pattern
                                                        } //  else LSC seq not found on =ve strand so ignore
                                                    } // no consenses seq found
                                                } // cov check
                                                left_seq.clear();
                                                last_leftSC = 0;
                                                consective_leftSC = 0;
                                            } // consecutive > 3 check
                                        }
                                    }
                                    left_seq.clear();
                                    last_leftSC = 0;
                                    consective_leftSC = 0;
                                    
                                    if (!problem_identified) { // for given RSC with found seq no LSC found, its single-sided-SC (chimera pattern)
                                        int start_rsc, end_rsc;

                                        if (last_rightSC < found_right - 10) { // SC seq fond on right side of contig{ // 16june2021
                                            query_frag = contig_seq.substr(last_rightSC, found_right - last_rightSC); // middle suspicious fragment
                                            first_frag = contig_seq.substr(0, last_rightSC - 1); // fragment before suspicious fragment
                                            
                                            if (found_right + 1 < contig_seq.size())
                                                sec_frag = contig_seq.substr(found_right + 1, contig_seq.size() - 1); // fragment after
                                            else
                                                sec_frag = "";
                                            start_rsc = last_rightSC, end_rsc = found_right;
                                            
                                        } else {//if (last_rightSC >= found_right) { // SC seq fond on left side of contig{ // commented 27/12
//                                            end_rsc = last_rightSC, start_rsc = found_right + right_cons_seq.size();
//                                            query_frag = contig_seq.substr(start_rsc, end_rsc - start_rsc); // middle suspicious fragment
//                                            problem_identified = true; // translocationcut from start to end and put before start
//                                            sc_patterns.push_back(start_rsc);
//                                            sc_patterns.push_back(end_rsc);
//                                            sc_patterns.push_back(found_right);
//                                            out_file << contig_name << "\t" << "translocation" << "\t" << start_rsc << "-" << end_rsc << "\t" << "forward" << "\t" << "insert at:" << found_right << endl;
//                                            right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
//                                            right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
//                                            mm = mm - consective_rightSC;
                                            goto out_right_next;
                                        }

                                        //******************************************************************     
                                        // come here only when last_rightSC < found_right 
                                        string target = first_frag + string(query_frag.size(), 'N') + sec_frag;
                                        string blast_hit = perform_blast(query_frag, target, "");
                                        vector<string> temp2;// = utils.split(blast_hit, "\t");
                                        utils.str_split(blast_hit, temp2, delimiter);
                                        blast_found = temp2[0];
                                        int target_start = atoi(temp2[1].c_str());
                                        int target_end = atoi(temp2[2].c_str());

                                        bool repeat = false;
                                        
                                        if (!sc_patterns.empty()) { // check miss-assembled regions
                                            
                                            for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                
                                                if (abs(target_start - *sc) < repeat_TH) {
                                                    
                                                    for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                        
                                                        if ((abs(target_end - *sc) < repeat_TH)) {
                                                            
                                                            repeat = true;
                                                            break;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        if (repeat) {
                                            //                                            sc_patterns.push_back(target_start);
                                            //                                            sc_patterns.push_back(target_end);
                                            goto out_right_next;
                                            
                                        } else if (blast_found != "*") { // hit found
                                            
                                            sc_patterns.push_back(target_start);
                                            sc_patterns.push_back(target_end);

                                            if (blast_found == "+") { // hit found on positive strand
                                                
                                                problem_identified = true; //selfChimera
                                                sc_patterns.push_back(start_rsc);
                                                sc_patterns.push_back(end_rsc);
                                                sc_patterns.push_back(found_right);
                                                right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                mm = mm - consective_rightSC;
                                                out_file << contig_name << "\t" << "selfChimera" << "\t" << start_rsc << "-" << end_rsc << "\t" << "forward" << "\t" << "Discard" << endl;
                                                goto out_right_next;
                                                
                                            } else { // RC of query seq
                                                query_frag = utils.Rcomplement(query_frag);
                                                //blast_found = perform_blast(query_frag, first_frag, sec_frag);
                                                problem_identified = true; //selfChimera
                                                sc_patterns.push_back(start_rsc);
                                                sc_patterns.push_back(end_rsc);
                                                sc_patterns.push_back(found_right);
                                                right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                mm = mm - consective_rightSC;
                                                out_file << contig_name << "\t" << "selfChimera" << "\t" << start_rsc << "-" << end_rsc << "\t" << "Reverse complement" << "\t" << "Discard" << endl;
                                                goto out_right_next;
                                            }
                                        } else { // No hit found so check any SC pattern initiate in emitted region
                                            // start = lsc, end= rsc
                                            float lsc_occurance = 0;

                                            bool right_found_in_left = false;

                                            for (it = left_SC_base_pos.begin(); it != left_SC_base_pos.end(); it++) {

                                                if (abs(*it - start_rsc) <= left_right_dist) {
                                                    right_found_in_left = true;
                                                    break;
                                                }
                                            }

                                           // it = find(left_SC_base_pos.begin(), left_SC_base_pos.end(), start_rsc);

                                            if(right_found_in_left) { //if (it != left_SC_base_pos.end()) {
                                                
                                                int pos = std::distance(left_SC_base_pos.begin(), it);
                                                vector <string> seq_lsc_occ;
                                                
                                                for (int l_found = pos; left_SC_base_pos[l_found] != 0; l_found++) {
                                                    seq_lsc_occ.push_back(left_SC_seq[l_found]);
                                                    lsc_occurance++;
                                                }
                                                if (lsc_occurance >= consective_SC_TH && start_rsc != previous_last_leftSC) { // lsc pattern found at the start of emitted sequence
                                                    
                                                    int pos_check;
                                                    
                                                    if (cov.at(start_rsc - 1) > 1) // to avoid zero cov at sc position
                                                        pos_check = start_rsc ;
                                                    else 
                                                        pos_check = start_rsc - 1;

                                                    if (float(lsc_occurance / cov.at(pos_check - 1)) * 100 >= SC_support_TH2) { // -1 because coverage starts from 0 index not one
                                                        string cons_lsc_occ = utils.consensus_seq_CAP3(seq_lsc_occ, toleft, PATH);//"left");
                                                        
                                                        if (cons_lsc_occ.size() >= cons_seq_size_TH) {

                                                            string cons_lsc_blast = perform_blast(cons_lsc_occ, contig_seq); //boost::regex_search(contig_seq, present, regexPattern);
                                                            
                                                            vector<string> temp3;// = utils.split(cons_lsc_blast, "\t");
                                                            utils.str_split(cons_lsc_blast, temp3, delimiter);
                                                            int cons_left_found = atoi(temp3[0].c_str());
                                                            string cons_left_strand = temp3[1];
                                                            
                                                            if (cons_left_found != -1) {
                                                                
                                                                if (cons_left_strand == "+") {
                                                                    
                                                                    cons_left_found = cons_left_found + cons_lsc_occ.length() - 1; // end pos of match to insert
                                                                    problem_identified = true; // *translocation*
                                                                    sc_patterns.push_back(start_rsc);
                                                                    sc_patterns.push_back(end_rsc);
                                                                    sc_patterns.push_back(cons_left_found);
                                                                    
                                                                    out_file << contig_name << "\t" << "translocation" << "\t" << start_rsc << "-" << end_rsc << "\t" << "forward" << "\t" << "insert at: " << cons_left_found << endl;
                                                                    // erase treated left sc, so that won't check them later
                                                                    left_st_ind = pos;
                                                                    left_en_ind = pos + lsc_occurance;
                                                                    left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                    left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                    right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                    right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                    mm = mm - consective_rightSC;
                                                                    goto out_right_next;
                                                                    
                                                                } else{// if (cons_left_strand == "-") {
                                                                    cons_lsc_occ = utils.Rcomplement(cons_lsc_occ); // take reverse complement of LSC seq 
                                                                    problem_identified = true; // *inversion and translocation*
                                                                    sc_patterns.push_back(start_rsc);
                                                                    sc_patterns.push_back(end_rsc);
                                                                    sc_patterns.push_back(cons_left_found);
                                                                    out_file << contig_name << "\t" << "inversion_translocation" << "\t" << start_rsc << "-" << end_rsc << "\t" << "Reverse complement" << "\t" << "insert at: " << cons_left_found << endl;
                                                                    left_st_ind = pos;
                                                                    left_en_ind = pos + lsc_occurance;
                                                                    left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                    left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                    right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                    right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                    mm = mm - consective_rightSC;
                                                                    goto out_right_next;
                                                                }
                                                            }
                                                        }// consenses check
                                                    } // cov check
                                                }
                                                lsc_occurance = 0;
                                            }
                                            if (!problem_identified) { // no LSC at start found so now check right corner of emitted sequence have any RSC pattern
                                                float rsc_occurance = 0;

                                                bool right_found_in_right = false;

                                                for (it = right_SC_base_pos.begin(); it != right_SC_base_pos.end(); it++) {

                                                    if (abs(*it - found_right) <= left_right_dist) {
                                                        right_found_in_right = true;
                                                        break;
                                                    }
                                                }

                                                //it = find(right_SC_base_pos.begin(), right_SC_base_pos.end(), found_right);

                                                if(right_found_in_right) { //if (it != right_SC_base_pos.end()) {
                                                    int pos = std::distance(right_SC_base_pos.begin(), it);
                                                    vector <string> seq_rsc_occ;

                                                    for (int r_found = pos; right_SC_base_pos[r_found] != 0; r_found++) {
                                                        seq_rsc_occ.push_back(right_SC_seq[r_found]);
                                                        rsc_occurance++;
                                                    }
                                                    if (rsc_occurance >= consective_SC_TH && rsc_occurance != previous_last_rightSC) { // rsc pattern found at the end of emitted sequence
                                                        int pos_check;
                                                       
                                                        if (cov.at(found_right - 1) > 1) // to avoid zero cov at sc position
                                                            pos_check = found_right ;
                                                        else 
                                                            pos_check = found_right - 1;
                                                        
                                                        if (float(rsc_occurance / cov.at(pos_check - 1)) * 100 >= SC_support_TH2) { // -1 because coverage starts from 0 index not one
                                                            string cons_rsc_occ = utils.consensus_seq_CAP3(seq_rsc_occ, !toleft, PATH); //"right");
                                                            
                                                            if (cons_rsc_occ.size() >= cons_seq_size_TH) {
                                                                
                                                                string cons_rsc_blast = perform_blast(cons_rsc_occ, contig_seq); // boost::regex_search(contig_seq, present, regexPattern);
                                                                vector<string> temp3;// = utils.split(cons_rsc_blast, "\t");
                                                                utils.str_split(cons_rsc_blast, temp3, delimiter);
                                                                int cons_right_found = atoi(temp3[0].c_str());
                                                                string cons_right_strand = temp3[1];
                                                                
                                                                if (cons_right_found != -1) {
                                                                    
                                                                    if (cons_right_strand == "+") {
                                                                        //int cons_right_found = present.position(a); //start pos of match to insert because forward // merge start and end
                                                                        problem_identified = true; // *translocation*
                                                                        sc_patterns.push_back(start_rsc);
                                                                        sc_patterns.push_back(end_rsc);
                                                                        sc_patterns.push_back(cons_right_found);
                                                                        out_file << contig_name << "\t" << "translocation" << "\t" << start_rsc << "-" << end_rsc << "\t" << "forward" << "\t" << "insert at: " << cons_right_found << endl;
                                                                        right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                        right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                        mm = mm - consective_rightSC;
                                                                        goto out_right_next;
                                                                        
                                                                    } else{// if (cons_right_strand == "-") {
                                                                        
                                                                        cons_rsc_occ = utils.Rcomplement(cons_rsc_occ); // take reverse complement of RSC seq 
                                                                        cons_right_found = cons_right_found + cons_rsc_occ.size(); // start pos of match to insert because RC
                                                                        problem_identified = true; // *inversion and translocation*
                                                                        sc_patterns.push_back(start_rsc);
                                                                        sc_patterns.push_back(end_rsc);
                                                                        sc_patterns.push_back(cons_right_found);
                                                                        out_file << contig_name << "\t" << "inversion_translocation" << "\t" << start_rsc << "-" << end_rsc << "\t" << "Reverse complement" << "\t" << "insert at: " << cons_right_found << endl;
                                                                        right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                        right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                        mm = mm - consective_rightSC;
                                                                        goto out_right_next;
                                                                    }
                                                                }
                                                            } // consesne check
                                                        } // cov check
                                                    }
                                                    rsc_occurance = 0;
                                                }
                                            }
                                            //added after meeting 3/12/2019
                                            if (!problem_identified) { // simply do the translocation of part containg RSC to the pos where it found, translocation better than false chimera
                                                first_frag = contig_seq.substr(0, last_rightSC);
                                                int start, end;
                                                sec_frag = contig_seq.substr(last_rightSC + 1, contig_seq.size() - 1); // fragment after

                                                if (first_frag.size() < sec_frag.size()) { // take smallest fragment of contig to check for chimera or local mis-assembly
                                                    query_frag = first_frag;
                                                    first_frag = "";
                                                    start = 0;
                                                    end = last_rightSC;
                                                    string blast_hit = perform_blast(query_frag, sec_frag, "");
                                                } else { // sec fragment is smaller
                                                    query_frag = sec_frag;
                                                    sec_frag = "";
                                                    start = last_rightSC;
                                                    end = contig_seq.length() - 1 ;
                                                    string blast_hit = perform_blast(query_frag, first_frag, "");
                                                }
                                                vector<string> temp2;// = utils.split(blast_hit, "\t");
                                                utils.str_split(blast_hit, temp2, delimiter);
                                                blast_found = temp2[0];
                                                int target_start = atoi(temp2[1].c_str());
                                                int target_end = atoi(temp2[2].c_str());
                                                bool repeat = false;
                                                
                                                if (!sc_patterns.empty()) { // check miss-assembled regions
                                                    
                                                    for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                        
                                                        if (abs(target_start - *sc) < repeat_TH) {
                                                            
                                                            for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                
                                                                if ((abs(target_end - *sc) < repeat_TH)) {
                                                                    
                                                                    repeat = true;
                                                                    break;
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                                if (repeat) {
                                                    //                                        sc_patterns.push_back(target_start);
                                                    //                                        sc_patterns.push_back(target_end);
                                                    goto out_right_next;
                                                }
                                                if (blast_found != "*") {
                                                    sc_patterns.push_back(target_start);
                                                    sc_patterns.push_back(target_end);
                                                    problem_identified = true; // selfChimera  start and end pos
                                                    sc_patterns.push_back(start);
                                                    sc_patterns.push_back(end);
                                                    sc_patterns.push_back(target_start);
                                                    sc_patterns.push_back(target_end);
                                                    
                                                    if (blast_found == "+") {
                                                        out_file << contig_name << "\t" << "selfChimera" << "\t" << start << "-" << end << "\t" << "Forward" << "\t" << "Discard" << endl;
                                                    } else {
                                                        query_frag = utils.Rcomplement(query_frag);
                                                        out_file << contig_name << "\t" << "selfChimera" << "\t" << start << "-" << end << "\t" << "Reverse complement" << "\t" << "Discard" << endl;
                                                    }
                                                    
                                                    problem_identified = true;
                                                    goto out_right_next;

                                                } else {
                                                    int start_rsc = last_rightSC, end_rsc = found_right;
                                                    query_frag = contig_seq.substr(start_rsc, end_rsc - start_rsc); // join two points and emit middle part // when last_leftSC > found_left
                                                    int sum = 0;
                                                    for (int aa = start_rsc - 110; aa <= start_rsc - 10; aa++) {
                                                        sum = sum + cov[aa];
                                                    }
                                                    int left_avg_cov = sum / 100;
                                                    sum = 0;

                                                    for (int aa = end_rsc + 10; aa <= end_rsc + 110; aa++) {
                                                        sum = sum + cov[aa];
                                                    }
                                                    int right_avg_cov = sum / 100;
                                                    sum = 0;

                                                    for (int aa = start_rsc; aa <= end_rsc; aa++) {
                                                        sum = sum + cov[aa];
                                                    }
                                                    int query_frag_avg = float(sum) / float(query_frag.size());
                                                    bool repeat = false;
                                                    
                                                    if (abs(end_rsc - start_rsc) >= ignore_short_seq && query_frag_avg >= query_frag_avg_TH) {
                                                        
                                                        if (!sc_patterns.empty()) { // check miss-assembled regions
                                                            
                                                            for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                
                                                                if (abs(start_rsc - *sc) < repeat_TH) {
                                                                    
                                                                    for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                        
                                                                        if ((abs(end_rsc - *sc) < repeat_TH)) {
                                                                            repeat = true;
                                                                            break;
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                        if (repeat) {
                                                            goto out_right_next;
                                                        }

                                                            //28oct21 avoid chimeras at the start or end positions of overlapped merged positions found in RI and merged by SCs
                                                                 float SC_TH_mc;

                                                        if (find(left_SC_base_pos.begin(), left_SC_base_pos.end(), start_rsc  - 1) != left_SC_base_pos.end()) { // if crisscross pattern
                                                            SC_TH_mc = SC_support_TH2;
                                                        } else // single sided pattern
                                                            SC_TH_mc = SC_support_TH;
                                                        // cout <<float(consective_rightSC / cov.at(pos_check - 1)) << endl;
                                                        if (float(consective_rightSC / cov.at(pos_check - 1)) * 100 >= SC_TH_mc) { // 10feb22
                                                           
                                                        if (merged_contigs_itr == merged_contigs.end()) {
                                                            
                                                            problem_identified = true; //merge start and end, trasns chimera 
                                                            sc_patterns.push_back(start_rsc);
                                                            sc_patterns.push_back(end_rsc);
                                                            out_file << contig_name << "\t" << "multigeneChimera" << "\t" << start_rsc << "-" << end_rsc << "\t" << "Forward" << "\t" << "New contig " << endl;
                                                            right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                            right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                            mm = mm - consective_rightSC;
                                                       
                                                        } else { //identified chimera positions must be at one side of merged points in the merged contigs

                                                            if ((start_rsc < start_overlap_merged && end_rsc < start_overlap_merged) || (start_rsc > end_overlap_merged && end_rsc > end_overlap_merged)) {
                                                                problem_identified = true; //merge start and end, trasns chimera 
                                                                sc_patterns.push_back(start_rsc);
                                                                sc_patterns.push_back(end_rsc);
                                                                out_file << contig_name << "\t" << "multigeneChimera" << "\t" << start_rsc << "-" << end_rsc << "\t" << "Forward" << "\t" << "New contig " << endl;
                                                                right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                mm = mm - consective_rightSC;

                                                            } // else go out
                                                        }
                                                        }
                                                            goto out_right_next;

                                                            } else if (query_frag.size() < ignore_short_seq) {
                                                        problem_identified = true; // merge start and end discard middle part
                                                        sc_patterns.push_back(start_rsc);
                                                        sc_patterns.push_back(end_rsc);
                                                        out_file << contig_name << "\t" << "unsupported_insertion" << "\t" << start_rsc << "-" << end_rsc<< "\t" << "Forward" << "\t" << "Discard" << endl;
                                                        right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                        right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                        mm = mm - consective_rightSC;
                                                        goto out_right_next;
                                                    } else {
                                                        sc_patterns.push_back(end_rsc);
                                                        sc_patterns.push_back(start_rsc);
                                                        out_file << contig_name << "\t" << "translocation" << "\t" << 0 << "-" << start << "\t" << "forwrd" << "\t" << "insert at: " << end << endl;
                                                        right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                        right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                        mm = mm - consective_rightSC;
                                                        problem_identified = true;
                                                        goto out_right_next;
                                                    }
                                                }
                                            }
                                            //**********************************************************************
                                            //}
                                        }
                                        //                                    if (!problem_identified) { // RSC found no LSC and no other pattern for selfChimera for omitted fragment
                                        //                                        int start_rsc = last_rightSC, end_rsc = found_right;
                                        //                                        query_frag = contig_seq.substr(start_rsc, end_rsc - start_rsc); // join two points and emit middle part // when last_leftSC > found_left
                                        //                                        int sum;
                                        //                                        for (int aa = start_rsc - 110; aa <= start_rsc - 10; aa++) {
                                        //                                            sum = sum + cov[aa];
                                        //                                        }
                                        //                                        int left_avg_cov = sum / 100;
                                        //                                        sum = 0;
                                        //
                                        //                                        for (int aa = end_rsc + 10; aa <= end_rsc + 110; aa++) {
                                        //                                            sum = sum + cov[aa];
                                        //                                        }
                                        //                                        int right_avg_cov = sum / 100;
                                        //                                        sum = 0;
                                        //
                                        //                                        for (int aa = start_rsc; aa <= end_rsc; aa++) {
                                        //                                            sum = sum + cov[aa];
                                        //                                        }
                                        //                                        int query_frag_avg = (float(sum) / float(query_frag.size()))*100;
                                        //                                        bool repeat = false;
                                        //                                        if (query_frag.size() >= 100 && query_frag_avg >= 25) {
                                        //                                            if (!sc_patterns.empty()) { // check miss-assembled regions
                                        //                                                for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                        //                                                    if (abs(start_rsc - *sc) < 100) {
                                        //                                                        for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                        //                                                            if ((abs(end_rsc - *sc) < 100)) {
                                        //                                                                repeat = true;
                                        //                                                                break;
                                        //                                                            }
                                        //                                                        }
                                        //                                                    }
                                        //                                                }
                                        //                                            }
                                        //                                            if (repeat == true) {
                                        //                                                goto out_right_next;
                                        //                                            }
                                        //                                            problem_identified = true; //merge start and end, trasns chimera 
                                        //                                            sc_patterns.push_back(start_rsc);
                                        //                                            sc_patterns.push_back(end_rsc);
                                        //                                            out_file << contig_name << "\t" << "multigeneChimera" << "\t" << start_rsc << "-" << end_rsc << "\t" << "Forward" << "\t" << "New contig " << endl;
                                        //                                            right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                        //                                            right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                        //                                            mm = mm - consective_rightSC;
                                        //                                            goto out_right_next;
                                        //                                        } else if (query_frag.size() < ignore_short_seq){  // 17 may 2021 bug fixed {
                                        //                                            problem_identified = true; // merge start and end discard middle part
                                        //                                            sc_patterns.push_back(start_rsc);
                                        //                                            sc_patterns.push_back(end_rsc);
                                        //                                            out_file << contig_name << "\t" << "unsupported_insertion" << "\t" << start_rsc << "-" << end_rsc << "\t" << "Forward" << "\t" << "Discard" << endl;
                                        //                                            right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                        //                                            right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                        //                                            mm = mm - consective_rightSC;
                                        //                                            goto out_right_next;
                                        //                                        }
                                    }
                                } else{// if (strand == "-") { // check for RC RSC
                                    //******************************************************************************************************************************************************                 
                                    found_right_end = found_right - atoi(temp[2].c_str()) ; // because st 20 end 10
                                    right_cons_seq = utils.Rcomplement(right_cons_seq);
                                    found_right = found_right;// + right_cons_seq.size(); // for right start pos of match should be same as start of LSC (start of read containing LSC)

                                    // check for towards pattern if RSc seq found somewhere

                                    for (int kk = 0; kk < left_SC_base_pos.size(); kk++) {
                                        //   cout<<right_SC_base_pos[kk]<<endl;
                                        if (left_SC_base_pos[kk] != 0) {
                                            last_leftSC = left_SC_base_pos[kk];
                                            consective_leftSC += 1;
                                            left_seq.push_back(left_SC_seq[kk]);

                                        } else {
                                            
                                            if (consective_leftSC >= consective_SC_TH && last_leftSC != previous_last_leftSC) {// && max_consec_rightSC < consective_rightSC) //max cluster supported by  >4 SC reads
                                                int pos_check;
                                                
                                                if (cov.at(last_leftSC - 1 ) > 1) // to avoid zero cov at sc position
                                                    pos_check = last_leftSC;
                                                else 
                                                    pos_check = last_leftSC - 1;
                                                
                                                if (float(consective_leftSC / cov.at(pos_check - 1)) * 100 >= SC_support_TH2) { // -1 because coverage starts from 0 index not one
                                                    left_st_ind = kk - consective_leftSC;
                                                    left_en_ind = kk - 1;
                                                    left_cons_seq = utils.consensus_seq_CAP3(left_seq, toleft, PATH);//"left");
                                                    
                                                    if (left_cons_seq.size() >= cons_seq_size_TH) {
                                                        previous_last_leftSC = last_leftSC;
                                                      //  left_cons_seq = utils.Rcomplement(left_cons_seq); // because RSC mapped in RC order so the LSC should too RC

                                                        string blast_hit_left = perform_blast(left_cons_seq, contig_seq); //boost::regex_search(contig_seq, present, regexPattern);
                                                        vector<string> temp2;// = utils.split(blast_hit_left, "\t");
                                                        utils.str_split(blast_hit_left, temp2, delimiter);
                                                        found_left = atoi(temp2[0].c_str());
                                                        found_left_end = atoi(temp2[2].c_str()); // 17jan - found_left;
                                                        string strand = temp2[1];

                                                        if (found_left != -1 && abs(found_left - last_leftSC) > left_right_dist && strand == "-") { //if (found_left != -1 && strand == "-") {

                                                            if ((abs(found_right - last_leftSC) <= left_right_dist) && (abs((found_left_end) - last_rightSC) <= left_right_dist)) { //towards pattern, seq match may not at exact LSC start 

                                                                int start_region = last_rightSC;
                                                                int end_region = found_right;
                                                                //*******************************************************************************************
                                                                if (end_region > start_region) {

                                                                    query_frag = contig_seq.substr(start_region, end_region - start_region); // middle suspicious fragment
                                                                    first_frag = contig_seq.substr(0, start_region - 1); // fragment before suspicious fragment

                                                                    if (end_region + 1 < contig_seq.size())
                                                                        sec_frag = contig_seq.substr(end_region + 1, contig_seq.size() - 1); // fragment after
                                                                    else
                                                                        sec_frag = "";
                                                             } else if (end_region < start_region) { // second pattern of RSC get match at first pattern of LSC 5feb22

                                                                    query_frag = contig_seq.substr(end_region, start_region - end_region); // middle suspicious fragment
                                                                    first_frag = contig_seq.substr(0, end_region - 1); // fragment before suspicious fragment

                                                                    if (start_region + 1 < contig_seq.size())
                                                                        sec_frag = contig_seq.substr(start_region + 1, contig_seq.size() - 1); // fragment after
                                                                    else
                                                                        sec_frag = "";

                                                                }
                                                                string target = first_frag + string(query_frag.size(), 'N') + sec_frag; // to obtain correct coordinates of hit

                                                                string blast_hit = perform_blast(query_frag, target, "");
                                                                vector<string> temp2; // = utils.split(blast_hit, "\t");
                                                                utils.str_split(blast_hit, temp2, delimiter);
                                                                blast_found = temp2[0];
                                                                int target_start = atoi(temp2[1].c_str());
                                                                int target_end = atoi(temp2[2].c_str());
                                                                bool repeat = false;

                                                                if (!sc_patterns.empty()) { // check miss-assembled regions

                                                                    for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {

                                                                        if (abs(target_start - *sc) < repeat_TH) {
                                                                            
                                                                            for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                                
                                                                                if ((abs(target_end - *sc) < repeat_TH)) {
                                                                                    repeat = true;
                                                                                    break;
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                                if (repeat) {
                                                                    //                                                                    sc_patterns.push_back(target_start);
                                                                    //                                                                    sc_patterns.push_back(target_end);
                                                                    goto out_right_next;
                                                                } else if (blast_found != "*") {
                                                                    sc_patterns.push_back(target_start);
                                                                    sc_patterns.push_back(target_end);
                                                                    
                                                                    if (blast_found == "+") {
                                                                        problem_identified = true; //selfChimera with inversion and translocation
                                                                        sc_patterns.push_back(start_region);
                                                                        sc_patterns.push_back(end_region);
                                                                    //27/12 edited because, it -- = + not possible to RC bboth side, SC won't match in thi case as RC RSC match to LSC Not RC of LSC and vice versa 
                                                                        //     out_file << contig_name << "\t" << "inversion_translocation&selfChimera" << "\t" << 0 << "-" << start_region << "\t" << "at:" << end_region << "\t" << end_region << "-" << contig_seq.length() - 1 << "\t" << "forward" << "\t" << "insert at:" << start_region << "\t" << start_region << "-" << end_region << "\t" << "discard" << endl;
                                                                        out_file << contig_name << "\t" << "selfChimera" << "\t" << start_region << "-" << end_region << "\t" << "Forward" << "\t" << "Discard" << endl;
                                                                        left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                        left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                        kk = kk - consective_leftSC;
                                                                        right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                        right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                        mm = mm - consective_rightSC;
                                                                        goto out_right_next;
                                                                        
                                                                    } else { // RC of query seq
                                                                        query_frag = utils.Rcomplement(query_frag);
                                                                        //blast_found = perform_blast(query_frag, first_frag, sec_frag);
                                                                        problem_identified = true; //selfChimera with inversion and translocation
                                                                        sc_patterns.push_back(start_region);
                                                                        sc_patterns.push_back(end_region);
                                                                       //out_file << contig_name << "\t" << "inversion_translocation&selfChimera" << "\t" << 0 << "-" << start_region << "\t" << "at:" << end_region << "\t" << end_region << "-" << contig_seq.length() - 1 << "\t" << "reverse complement" << "\t" << "insert at:" << start_region << "\t" << start_region << "-" << end_region << "\t" << "discard" << endl;
                                                                        out_file << contig_name << "\t" << "selfChimera" << "\t" << start_region << "-" << end_region << "\t" << "reverse" << "\t" << "Discard" << endl;
                                                                        left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                        left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                        kk = kk - consective_leftSC;
                                                                        right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                        right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                        mm = mm - consective_rightSC;
                                                                        goto out_right_next;
                                                                    }
                                                                } //8/1/2020
                                                                else {
                                                                    problem_identified = true; // unsupported_insertion, join two points and discard middle part
                                                                    sc_patterns.push_back(found_right);
                                                                    sc_patterns.push_back(last_rightSC);
                                                                    out_file << contig_name << "\t" << "inversion" << "\t"  << last_rightSC << "-" << found_right   << "\t" << "reverse complement" << "\t" << "invert at its own place" << endl;
                                                                    right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                    right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                    mm = mm - consective_rightSC;
                                                                    goto out_right_next;
                                                                }
                                                            }
                                                        } // else LSC seq not found
                                                    }//consenses check
                                                } // cov check
                                            } // no LSC found for towards pattern at all 
                                            left_seq.clear();
                                            last_leftSC = 0;
                                            consective_leftSC = 0;
                                        }
                                    }
                                    left_seq.clear();
                                    last_leftSC = 0;
                                    consective_leftSC = 0;
                                    
                                    if (!problem_identified) { // for given RC of RSC with found seq no LSC found and RSC found, //its single-sided-SC (chimera pattern)
                                        int start_rsc, end_rsc;
                                      
                                        if (last_rightSC < found_right -10 ) { // 16 june 2021
                                            query_frag = contig_seq.substr(last_rightSC, found_right - last_rightSC); // middle suspicious fragment
                                            first_frag = contig_seq.substr(0, last_rightSC - 1); // fragment before suspicious fragment
                                        
                                            if (found_right + 1 < contig_seq.size()) // hit found end of contig
                                                sec_frag = contig_seq.substr(found_right + 1, contig_seq.size() - 1); // fragment after
                                            else
                                                sec_frag = "";
                                            start_rsc = last_rightSC, end_rsc = found_right;

                                        } else if (last_rightSC > found_right + found_right_end) { // SC seq fond on left side of contig{ // edited 28/12/19
                                            end_rsc = last_rightSC - 1;   // 6jan2021
                                            start_rsc = found_right + found_right_end ;
                                            query_frag = contig_seq.substr(start_rsc, end_rsc - start_rsc);
                                            problem_identified = true; // Inversion from start_rsc to end_rsc 
                                            sc_patterns.push_back(start_rsc);
                                            sc_patterns.push_back(end_rsc);
                                            out_file << contig_name << "\t" << "inversion" << "\t" << start_rsc << "-" << end_rsc << "\t" << "reverse complement" << "\t" << "invert at its own place" << endl;
                                            right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                            right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                            mm = mm - consective_rightSC;
                                            goto out_right_next;
                                        }
                                        else
                                            goto out_right_next;

                                        //******************************************************************     
                                        // come here only when last_rightSC < found_right
                                        string target = first_frag + string(query_frag.size(), 'N') + sec_frag;
                                        string blast_hit = perform_blast(query_frag, target, "");
                                        vector<string> temp2;// = utils.split(blast_hit, "\t");
                                        utils.str_split(blast_hit, temp2, delimiter);
                                        blast_found = temp2[0];
                                        int target_start = atoi(temp2[1].c_str());
                                        int target_end = atoi(temp2[2].c_str());
                                        bool repeat = false;
                                        
                                        if (!sc_patterns.empty()) { // check miss-assembled regions
                                            
                                            for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                
                                                if (abs(target_start - *sc) < repeat_TH) {
                                                    
                                                    for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                        
                                                        if ((abs(target_end - *sc) < repeat_TH)) {
                                                            repeat = true;
                                                            break;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        if (repeat) {
                                            //                                            sc_patterns.push_back(target_start);
                                            //                                            sc_patterns.push_back(target_end);
                                            goto out_right_next;
                                        } else if (blast_found != "*") {
                                            sc_patterns.push_back(target_start);
                                            sc_patterns.push_back(target_end);
                                            
                                            if (blast_found == "+") {
                                                problem_identified = true; //selfChimera with inversion and translocation
                                                sc_patterns.push_back(start_rsc);
                                                sc_patterns.push_back(end_rsc);
                                                sc_patterns.push_back(found_right);
                                                //edited 27/12
                                                out_file << contig_name << "\t" << "selfChimera" << "\t" << start_rsc << "-" << end_rsc << "\t" << "Forward" << "\t" << "Discard" << endl;
                                               // out_file << contig_name << "\t" << "inversion_translocation&selfChimera" << "\t" << 0 << "-" << start_rsc << "\t" << "at:" << end_rsc << "\t" << end_rsc << "-" << contig_seq.length() - 1 << "\t" << "reverse complement" << "\t" << "insert at:" << start_rsc << "\t" << start_rsc << "-" << end_rsc << "\t" << "discard" << endl;
                                                right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                mm = mm - consective_rightSC;
                                                goto out_right_next;

                                            } else { // RC of query seq
                                                query_frag = utils.Rcomplement(query_frag);
                                                //blast_found = perform_blast(query_frag, first_frag, sec_frag);
                                                problem_identified = true; //selfChimera with inversion and translocation
                                                sc_patterns.push_back(start_rsc);
                                                sc_patterns.push_back(end_rsc);
                                                sc_patterns.push_back(found_right);
                                               //edited 27/12
                                                out_file << contig_name << "\t" << "selfChimera" << "\t" << start_rsc << "-" << end_rsc << "\t" << "Forward" << "\t" << "Discard" << endl;
                                                // out_file << contig_name << "\t" << "inversion_translocation&selfChimera" << "\t" << 0 << "-" << start_rsc << "\t" << "at:" << end_rsc << "\t" << end_rsc << "-" << contig_seq.length() - 1 << "\t" << "reverse complement" << "\t" << "insert at:" << start_rsc << "\t" << start_rsc << "-" << end_rsc << "\t" << "discard" << endl;
                                                right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                mm = mm - consective_rightSC;
                                                goto out_right_next;
                                            }
                                        } else { // No hit found so check any SC pattern initiate in emitted region          
//                                            query_frag = contig_seq.substr(0, last_rightSC); //insert at found_right + cons_seq // when last_rightSC < found_right and RC of RSC found 
//                                            problem_identified = true; // simple inversion and translocationof query fragment, insert
//                                            sc_patterns.push_back(last_rightSC);
//                                            sc_patterns.push_back(found_right);
//                                            out_file << contig_name << "\t" << "inversion_translocation" << "\t" << 0 << "-" << last_rightSC << "\t" << "reverse complement" << "\t" << "inversion and translocationat:" << found_right << endl;
//                                            right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
//                                            right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
//                                            mm = mm - consective_rightSC;
                                              //int start_rsc = last_rightSC, end_rsc = found_right;
                                                 //   query_frag = contig_seq.substr(start_rsc, end_rsc - start_rsc); // join two points and emit middle part // when last_leftSC > found_left
                                                    int sum = 0;
                                                    for (int aa = start_rsc - 110; aa <= start_rsc - 10; aa++) {
                                                        sum = sum + cov[aa];
                                                    }
                                                    int left_avg_cov = sum / 100;
                                                    sum = 0;

                                                    for (int aa = end_rsc + 10; aa <= end_rsc + 110; aa++) {
                                                        sum = sum + cov[aa];
                                                    }
                                                    int right_avg_cov = sum / 100;
                                                    sum = 0;

                                                    for (int aa = start_rsc; aa <= end_rsc; aa++) {
                                                        sum = sum + cov[aa];
                                                    }
                                                    int query_frag_avg = float(sum) / float(query_frag.size());
                                                    bool repeat = false;
                                                    
                                                    if (abs(end_rsc - start_rsc) >= ignore_short_seq && query_frag_avg >= query_frag_avg_TH) {
                                                        
                                                        if (!sc_patterns.empty()) { // check miss-assembled regions
                                                            
                                                            for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                
                                                                if (abs(start_rsc - *sc) < repeat_TH) {
                                                                    
                                                                    for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                        
                                                                        if ((abs(end_rsc - *sc) < repeat_TH)) {
                                                                            repeat = true;
                                                                    break;
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                                if (repeat) {
                                                    goto out_right_next;
                                                    }

                                                    //28oct21 avoid chimeras at the start or end positions of overlapped merged positions found in RI and merged by SCs
                                                float SC_TH_mc;

                                                if (find(left_SC_base_pos.begin(), left_SC_base_pos.end(), start_rsc - 1) != left_SC_base_pos.end()) { // if crisscross pattern
                                                    SC_TH_mc = SC_support_TH2;
                                                } else // single sided pattern
                                                    SC_TH_mc = SC_support_TH;
 //cout <<float(consective_rightSC / cov.at(pos_check - 1)) << endl;
                                                if (float(consective_rightSC / cov.at(pos_check - 1)) * 100 >= SC_TH_mc) { // 10feb22

                                                if (merged_contigs_itr == merged_contigs.end()) {
                                                    problem_identified = true; //merge start and end, trasns chimera 
                                                    sc_patterns.push_back(start_rsc);
                                                    sc_patterns.push_back(end_rsc);
                                                    out_file << contig_name << "\t" << "multigeneChimera" << "\t" << start_rsc << "-" << end_rsc << "\t" << "Forward" << "\t" << "New contig " << endl;
                                                    right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                    right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                    mm = mm - consective_rightSC;
                                                    
                                                } else { //identified chimera positions must be at one side of merged points in the merged contigs
                                                    
                                                    if ((start_rsc < start_overlap_merged && end_rsc < start_overlap_merged) || (start_rsc > end_overlap_merged && end_rsc > end_overlap_merged)) {

                                                        problem_identified = true; //merge start and end, trasns chimera 
                                                        sc_patterns.push_back(start_rsc);
                                                        sc_patterns.push_back(end_rsc);
                                                        out_file << contig_name << "\t" << "multigeneChimera" << "\t" << start_rsc << "-" << end_rsc << "\t" << "Forward" << "\t" << "New contig " << endl;
                                                        right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                        right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                        mm = mm - consective_rightSC;

                                                    } // else go out
                                                }
                                                }   
                                                        
                                                goto out_right_next;

                                            } else if (query_frag.size() < ignore_short_seq) {
                                                problem_identified = true; // merge start and end discard middle part
                                                sc_patterns.push_back(start_rsc);
                                                sc_patterns.push_back(end_rsc);
                                                out_file << contig_name << "\t" << "unsupported_insertion" << "\t" << start_rsc << "-" << end_rsc << "\t" << "Forward" << "\t" << "Discard" << endl;
                                                right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                mm = mm - consective_rightSC;
                                                goto out_right_next;
                                            } else {
                                                //                                                        sc_patterns.push_back(end_rsc);
                                                //                                                        sc_patterns.push_back(start_rsc);
                                                //                                                        out_file << contig_name << "\t" << "translocation" << "\t" << 0 << "-" << start_rsc << "\t" << "forwrd" << "\t" << "insert at: " << end_rsc << endl;
                                                //                                                        right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                //                                                        right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                //                                                        mm = mm - consective_rightSC;
                                                problem_identified = true;
                                                goto out_right_next;
                                            }
                                           // goto out_right_next;
                                        } // +ve strand
                                    }
                                }
                            } else { // RSC seq not found so check for crossed pattern
                                int start_lsc, end_lsc;
                                
                                for (int kk = 0; kk < left_SC_base_pos.size(); kk++) {
                                    //   cout<<right_SC_base_pos[kk]<<endl;
                                    if (left_SC_base_pos[kk] != 0) {
                                        last_leftSC = left_SC_base_pos[kk];
                                        consective_leftSC += 1;
                                        left_seq.push_back(left_SC_seq[kk]);

                                    } else {
                                        
                                        if (consective_leftSC > consective_SC_TH && last_leftSC != previous_last_leftSC)// && max_consec_rightSC < consective_rightSC) //max cluster supported by  >4 SC reads
                                        {
                                            left_st_ind = kk - consective_leftSC;
                                            left_en_ind = kk - 1;
                                            left_cons_seq = utils.consensus_seq_CAP3(left_seq, toleft, PATH); //"left");
                                            
                                            if (left_cons_seq.size() >= cons_seq_size_TH) {
                                                // to cover missing seq like TRINITY_DN20900_c0_g2
                                                if ((abs(last_leftSC - last_rightSC) <= left_right_dist) || (last_leftSC - last_rightSC >= 0 && last_leftSC < last_rightSC + right_cons_seq.length())) {
                                                    previous_last_leftSC = last_leftSC;

                                                    string found_left_blast = perform_blast(left_cons_seq, contig_seq); //boost::regex_search(contig_seq, present, regexPattern);
                                                    vector<string> temp2;// = utils.split(found_left_blast, "\t");
                                                    utils.str_split(found_left_blast, temp2, delimiter);
                                                    found_left = atoi(temp2[0].c_str());
                                                    found_left_end =  atoi(temp2[2].c_str()); //17 jan - found_left; 
                                                    string found_left_strand = temp2[1];
                                                    
                                                    if (found_left != -1 && abs(found_left_end - last_leftSC) > left_right_dist) { // LSC seq found in contig -> merge two points 
                                                        
                                                        if (found_left_strand == "+") {
                                                            //found_left = present.position(a);
                                                            //extract middle part
                                                            // do blast -> self both for and RC
                                                            // else find pos < lsc pos -> check any LSC at start pos -> forword?  translocation-> RC? inversion and translocation(rsc ALREADY CEHCKED)
                                                            // if pos > lsc simple translocationof fragment between pos -lsc insert at end of match found  
                                                            if (last_leftSC > found_left + 10) { // 16 june 2021 to avaoid processing of immediate hits
                                                                start_lsc = found_left_end ; // +1
                                                                end_lsc = last_leftSC - 1; // 6 jan 2021
                                                                query_frag = contig_seq.substr(start_lsc, end_lsc - start_lsc); // middle suspicious fragment
                                                                first_frag = contig_seq.substr(0, start_lsc - 1); // fragment before suspicious fragment
                                                                
                                                                if (end_lsc + 1 < contig_seq.size())
                                                                    sec_frag = contig_seq.substr(end_lsc + 1, contig_seq.size() - 1); // fragment after
                                                                else
                                                                    sec_frag = "";
                                                                
                                                            } else{// if (last_leftSC <= found_left) { // SC seq fond on left side of contig{
                                                                // edited 27/12 because 
//                                                                end_lsc = found_left, start_lsc = last_leftSC;
//                                                                query_frag = contig_seq.substr(start_lsc, end_lsc - start_lsc); // middle suspicious fragment
//                                                                problem_identified = true; // translocationcut from start to end and put at found_left + cons_seq.size();
//                                                                sc_patterns.push_back(start_lsc);
//                                                                sc_patterns.push_back(end_lsc);
//                                                                sc_patterns.push_back(found_left);
//                                                                out_file << contig_name << "\t" << "translocation" << "\t" << start_lsc << "-" << end_lsc << "\t" << "forward" << "\t" << "insert at:" << found_left + left_cons_seq.size() - 1 << endl;
//                                                                left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
//                                                                left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
//                                                                kk = kk - consective_leftSC;
//                                                                right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
//                                                                right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
//                                                                mm = mm - consective_rightSC;
                                                                  goto out_right_next;
                                                            }

                                                            //come here only when last_leftSC  > found_left
                                                            //start******************************************************************
                                                            string target = first_frag + string(query_frag.size(), 'N') + sec_frag;
                                                            string blast_hit = perform_blast(query_frag, target, "");
                                                            vector<string> temp2;// = utils.split(blast_hit, "\t");
                                                            utils.str_split(blast_hit, temp2, delimiter);
                                                            blast_found = temp2[0];
                                                            int target_start = atoi(temp2[1].c_str());
                                                            int target_end = atoi(temp2[2].c_str());
                                                            bool repeat = false;
                                                            
                                                            if (!sc_patterns.empty()) { // check miss-assembled regions
                                                                
                                                                for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                    
                                                                    if (abs(target_start - *sc) < repeat_TH) {
                                                                        
                                                                        for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                            
                                                                            if ((abs(target_end - *sc) < repeat_TH)) {
                                                                                repeat = true;
                                                                                break;
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                            if (repeat) {
                                                                //                                                                sc_patterns.push_back(target_start);
                                                                //                                                                sc_patterns.push_back(target_end);
                                                                goto out_right_next;
                                                            } else if (blast_found != "*") {
                                                                sc_patterns.push_back(target_start);
                                                                sc_patterns.push_back(target_end);
                                                                
                                                                if (blast_found == "+") {
                                                                    problem_identified = true; //selfChimera -> merge first and sec fragment
                                                                    sc_patterns.push_back(start_lsc);
                                                                    sc_patterns.push_back(end_lsc);
                                                                    out_file << contig_name << "\t" << "selfChimera" << "\t" << start_lsc << "-" << end_lsc << "\t" << "forward" << "\t" << "Discard" << endl;
                                                                    left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                    left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                    kk = kk - consective_leftSC;
                                                                    right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                    right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                    mm = mm - consective_rightSC;
                                                                    goto out_right_next;
                                                                    
                                                                } else { // RC of query seq
                                                                    query_frag = utils.Rcomplement(query_frag);
                                                                    //blast_found = perform_blast(query_frag, first_frag, sec_frag);
                                                                    problem_identified = true; //selfChimera -> merge first and sec fragment
                                                                    sc_patterns.push_back(start_lsc);
                                                                    sc_patterns.push_back(end_lsc);
                                                                    out_file << contig_name << "\t" << "selfChimera" << "\t" << start_lsc << "-" << end_lsc << "\t" << "Reverse complement" << "\t" << "Discard" << endl;
                                                                    left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                    left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                    kk = kk - consective_leftSC;
                                                                    right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                    right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                    mm = mm - consective_rightSC;
                                                                    goto out_right_next;
                                                                }
                                                            } else { // check any SC pattern initiate from emitted region (query fragment)
                                                                // left_found = lsc, right_found = rsc
                                                                float lsc_occurance = 0;
                                                               // it = find(left_SC_base_pos.begin(), left_SC_base_pos.end(), start_lsc);

                                                                bool left_found_in_left = false;

                                                                for (it = left_SC_base_pos.begin(); it != left_SC_base_pos.end(); it++) {

                                                                    if (abs(*it - start_lsc) <= left_right_dist) {
                                                                        left_found_in_left = true;
                                                                        break;
                                                                    }
                                                                }
                                                                if (left_found_in_left) {//if (it != left_SC_base_pos.end()) {
                                                                    int pos = std::distance(left_SC_base_pos.begin(), it);
                                                                    vector <string> seq_lsc_occ;

                                                                    for (int l_found = pos; left_SC_base_pos[l_found] != 0; l_found++) {
                                                                        seq_lsc_occ.push_back(left_SC_seq[l_found]);
                                                                        lsc_occurance++;
                                                                    }
                                                                    if (lsc_occurance >= consective_SC_TH && lsc_occurance != previous_last_leftSC) { // lsc pattern found at the start of emitted sequence
                                                                        int pos_check;
                                                                        
                                                                        if (cov.at(start_lsc ) > 1) // to avoid zero cov at sc position -2 for RSC none for LSC
                                                                            pos_check = start_lsc;
                                                                        else 
                                                                            pos_check = start_lsc - 1;
                                                                        
                                                                        if (float(lsc_occurance / cov.at(pos_check - 1)) * 100 >= SC_support_TH2) { // -1 because coverage starts from 0 index not one
                                                                            string cons_lsc_occ = utils.consensus_seq_CAP3(seq_lsc_occ, toleft, PATH); //"left");
                                                                            
                                                                            if (cons_lsc_occ.size() >= cons_seq_size_TH) {
                                                                                //
                                                                                string cons_left_blast = perform_blast(cons_lsc_occ, contig_seq); //boost::regex_search(contig_seq, present, regexPattern);
                                                                                vector<string> temp3; //= utils.split(cons_left_blast, "\t");
                                                                                utils.str_split(cons_left_blast, temp3, delimiter);
                                                                                int cons_left_found = atoi(temp3[0].c_str());
                                                                                int cons_left_end = atoi(temp3[2].c_str()) - cons_left_found;
                                                                                string con_left_strand = temp3[1];
                                                                                
                                                                                if (cons_left_found != -1) {
                                                                                    
                                                                                    if (con_left_strand == "+") {
                                                                                        
                                                                                        cons_left_found = cons_left_found + cons_left_end ; // end pos of match to insert
                                                                                        left_st_ind = pos;
                                                                                        left_en_ind = pos + lsc_occurance;
                                                                                        left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                                        left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                                        kk = kk - consective_leftSC;
                                                                                        problem_identified = true; // *translocation* insert query after cons_left_found and merge first/sec frags
                                                                                        sc_patterns.push_back(start_lsc);
                                                                                        sc_patterns.push_back(end_lsc);
                                                                                        sc_patterns.push_back(cons_left_found);
                                                                                        out_file << contig_name << "\t" << "translocation" << "\t" << start_lsc << "-" << end_lsc << "\t" << "forward" << "\t" << "insert at:" << cons_left_found << endl;
                                                                                        right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                                        right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                                        mm = mm - consective_rightSC;
                                                                                        goto out_right_next;
                                                                                        
                                                                                    } else if (con_left_strand == "-") {
                                                                                        {
                                                                                            cons_lsc_occ = utils.Rcomplement(cons_lsc_occ); // take reverse complement of LSC seq 
                                                                                            //
                                                                                            //                                                                                    int cons_left_found = perform_blast(cons_lsc_occ, contig_seq); //boost::regex_search(contig_seq, present, regexPattern);
                                                                                            //                                                                                    if (cons_left_found != -1) {
                                                                                            //int cons_left_found = present.position(a); // start pos of match to insert because RC
                                                                                            problem_identified = true; // *inversion and translocation* insert query after cons_left_found and merge first/sec frags
                                                                                            sc_patterns.push_back(start_lsc);
                                                                                            sc_patterns.push_back(end_lsc);
                                                                                            sc_patterns.push_back(cons_left_found);
                                                                                            out_file << contig_name << "\t" << "inversion_translocation" << "\t" << start_lsc << "-" << end_lsc << "\t" << "Reverse complement" << "\t" << "insert at:" << cons_left_found << endl; 
                                                                                           // out_file << contig_name << "\t" << "inversion_translocation" << "\t" << 0 << "-" << start_lsc << "\t" << "insert at:" << end_lsc << "\t" << start_lsc << "-" << contig_seq.size() - 1 << "\t" << "reverse complement" << "\t" << "insert at:" << cons_left_found << endl;
                                                                                            left_st_ind = pos;
                                                                                            left_en_ind = pos + lsc_occurance;
                                                                                            left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                                            left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                                            kk = kk - consective_leftSC;
                                                                                            right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                                            right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                                            mm = mm - consective_rightSC;
                                                                                            goto out_right_next;
                                                                                        }
                                                                                    }
                                                                                } // cons check
                                                                            } // cov check
                                                                        }
                                                                        lsc_occurance = 0;
                                                                    }
                                                                    // if problem not identified yet don't check right as its already been checked above
                                                                }
                                                                //}
                                                                //end ******************************************************************************

                                                                if (!problem_identified) { // check if query fragment meet criteria of multigeneChimera because LSC seq found somewhere in contig while checking for crossed pattern
                                                                    query_frag = contig_seq.substr(start_lsc, end_lsc - start_lsc); // join two points and emit middle part // when last_leftSC > found_left
                                                                    int sum = 0;
//                                                                    for (int aa = start_lsc - 110; aa <= start_lsc - 10; aa++) {
//                                                                        sum = sum + cov[aa];
//                                                                    }
//                                                                    int left_avg_cov = sum / 100;
//                                                                    sum = 0;
//
//                                                                    for (int aa = end_lsc + 10; aa <= end_lsc + 110; aa++) {
//                                                                        sum = sum + cov[aa];
//                                                                    }
//                                                                    int right_avg_cov = sum / 100;
//                                                                    sum = 0;

                                                                    for (int aa = start_lsc; aa <= end_lsc; aa++) {
                                                                        sum = sum + cov[aa];
                                                                    }
                                                                    int query_frag_avg = float(sum) / float(query_frag.size());

                                                                    if (abs(end_lsc - start_lsc) >= ignore_short_seq && query_frag_avg >= query_frag_avg_TH) {
                                                                        
                                                                        if (!sc_patterns.empty()) { // check miss-assembled regions
                                                                            
                                                                            for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                                
                                                                                if (abs(start_lsc - *sc) < repeat_TH) {
                                                                                    
                                                                                    for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                                        
                                                                                        if ((abs(end_lsc - *sc) < repeat_TH)) {
                                                                                            repeat = true;
                                                                                            break;
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                        if (repeat) {
                                                                            goto out_right_next;
                                                                        }
                                                                                                                              //28oct21 avoid chimeras at the start or end positions of overlapped merged positions found in RI and merged by SCs
                                                                            //identified chimera positions must be at one side of merged points in the merged contigs
                                                                        float SC_TH_mc;
                                                                          
                                                                        if (find(right_SC_base_pos.begin(), right_SC_base_pos.end(), last_leftSC - 1) != right_SC_base_pos.end()) { // if crisscross pattern
                                                                            SC_TH_mc = SC_support_TH2;
                                                                        } else // single sided pattern
                                                                            SC_TH_mc = SC_support_TH;
 //cout <<float(consective_leftSC / cov.at(pos_check - 1)) << endl;
                                                                        if (float(consective_leftSC / cov.at(pos_check - 1)) * 100 >= SC_TH_mc) { // 10feb22

                                                                        if (merged_contigs_itr == merged_contigs.end()) {
                                                                            problem_identified = true;
                                                                            left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                            left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                            kk = kk - consective_leftSC;
                                                                            sc_patterns.push_back(start_lsc);
                                                                            sc_patterns.push_back(end_lsc);
                                                                            out_file << contig_name << "\t" << "multigeneChimera" << "\t" << start_lsc << "-" << end_lsc << "\t" << "forward" << "\t" << "New contig" << endl;
                                                                            right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                            right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                            mm = mm - consective_rightSC;
                                                                        } else {

                                                                            if ((start_lsc < start_overlap_merged && end_lsc < start_overlap_merged) || (start_lsc > end_overlap_merged && end_lsc > end_overlap_merged)) {
                                                                                problem_identified = true;
                                                                                left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                                left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                                kk = kk - consective_leftSC;
                                                                                sc_patterns.push_back(start_lsc);
                                                                                sc_patterns.push_back(end_lsc);
                                                                                out_file << contig_name << "\t" << "multigeneChimera" << "\t" << start_lsc << "-" << end_lsc << "\t" << "forward" << "\t" << "New contig" << endl;
                                                                                right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                                right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                                mm = mm - consective_rightSC;
                                                                            } // else go out

                                                                        }
                                                                        }
                                                                        goto out_right_next;
                                                                        
                                                                    } else if (query_frag.size() < ignore_short_seq){  // 17 may 2021 bug fixed{
                                                                        problem_identified = true; 
                                                                        sc_patterns.push_back(start_lsc);
                                                                        sc_patterns.push_back(end_lsc);
                                                                        out_file << contig_name << "\t" << "unsupported_insertion" << "\t" << start_lsc << "-" << end_lsc << "\t" << "forward" << "\t" << "Discard" << endl;
                                                                        left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                        left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                        kk = kk - consective_leftSC;
                                                                        right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                        right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                        mm = mm - consective_rightSC;
                                                                        goto out_right_next;
                                                                    }
                                                                }

                                                            }
                                                        } else{// if (found_left_strand == "-") { // for LSC seq not found in contig search RC of LSC
                                                            left_cons_seq = utils.Rcomplement(left_cons_seq);
                                                            //start*****************************************************************************
                                                            /* if(match found)
                                                             * simple inversion and translocationstart of match
                                                             * if pos > extract lsc till pos and add at pos+cons_seq_length ( invert the main part at its own position) -> just inversion
                                                             * if pos < extract cons_Seq_len till end of cotig and add at pos - > inversion and translocation(leave middle part as it is)
                                                             * else do following thing 
                                                             */

                                                            if (last_leftSC > found_left + 10) { // 16 june 2021
                                                                start_lsc = found_left + left_cons_seq.size() + 1;
                                                                end_lsc = last_leftSC;
                                                                query_frag = contig_seq.substr(start_lsc, end_lsc - start_lsc); // middle suspicious fragment
                                                                first_frag = contig_seq.substr(0, start_lsc - 1); // fragment before suspicious fragment // fro blast check purpose
                                                                
                                                                if (end_lsc + 1 < contig_seq.size())
                                                                    sec_frag = contig_seq.substr(end_lsc + 1, contig_seq.size() - 1 ) ; // fragment after
                                                                else
                                                                    sec_frag = "";

                                                            } else { //if (last_leftSC < found_left) { // SC seq fond on right side of contig{
                                                                
                                                                end_lsc = found_left, start_lsc = last_leftSC;
                                                                query_frag = contig_seq.substr(start_lsc, end_lsc - start_lsc); // middle suspicious fragment
                                                                problem_identified = true; //  cut from start to end and invert at its own place
                                                                sc_patterns.push_back(start_lsc);
                                                                sc_patterns.push_back(end_lsc);
                                                                out_file << contig_name << "\t" << "inversion" << "\t" << start_lsc << "-" << end_lsc << "\t" << "Reverse complement" << "\t" << "invert at its own place" << endl;

                                                                left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                kk = kk - consective_leftSC;
                                                                right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                mm = mm - consective_rightSC;
                                                                goto out_right_next;
                                                            }
                                                            string target = first_frag + string(query_frag.size(), 'N') + sec_frag;
                                                            string blast_hit = perform_blast(query_frag, target, "");
                                                            vector<string> temp2;// = utils.split(blast_hit, "\t");
                                                            utils.str_split(blast_hit, temp2, delimiter);
                                                            blast_found = temp2[0];
                                                            int target_start = atoi(temp2[1].c_str());
                                                            int target_end = atoi(temp2[2].c_str());
                                                            bool repeat = false;
                                                            
                                                            if (!sc_patterns.empty()) { // check miss-assembled regions
                                                                
                                                                for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                    
                                                                    if (abs(target_start - *sc) < repeat_TH) {
                                                                        
                                                                        for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                            
                                                                            if ((abs(target_end - *sc) < repeat_TH)) {
                                                                                repeat = true;
                                                                                break;
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                            if (repeat) {
                                                                //                                                                sc_patterns.push_back(target_start);
                                                                //                                                                sc_patterns.push_back(target_end);
                                                                goto out_right_next;
                                                            } else if (blast_found != "*") {
                                                                
                                                                sc_patterns.push_back(target_start);
                                                                sc_patterns.push_back(target_end);
                                                                
                                                                if (blast_found == "+") {
                                                                    problem_identified = true; //selfChimera with inversion and translocationof left frag at start of match
                                                                    sc_patterns.push_back(start_lsc);
                                                                    sc_patterns.push_back(end_lsc);
                                                                    sc_patterns.push_back(found_left);
                                                                    
                                                                    out_file << contig_name << "\t" << "inversion_translocation_selfChimera" << "\t" << end_lsc << "-" << contig_seq.size() << "\t" << "at:" << found_left << "\t" << start_lsc << "-" << end_lsc << "\t" << "forward" << "\t" << "discard" << endl;
                                                                     //out_file << contig_name << "\t" << "selfChimera" << "\t" << start_lsc << "-" << end_lsc << "\t" << "forward" << "\t" << "Discard" << endl;
                                                                    left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                    left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                    kk = kk - consective_leftSC;
                                                                    right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                    right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                    mm = mm - consective_rightSC;
                                                                    goto out_right_next;
                                                                    
                                                                } else { // RC of query seq
                                                                    query_frag = utils.Rcomplement(query_frag);
                                                                    //blast_found = perform_blast(query_frag, first_frag, sec_frag);
                                                                    problem_identified = true; //selfChimera with inversion and translocationof left frag at start of match
                                                                    sc_patterns.push_back(start_lsc);
                                                                    sc_patterns.push_back(end_lsc);
                                                                    sc_patterns.push_back(found_left);
                                                                    out_file << contig_name << "\t" << "inversion_translocation_selfChimera" << "\t" << end_lsc << "-" << contig_seq.size() << "\t" << "at:" << found_left << " \t" << start_lsc << "-" << end_lsc << "\t" << "Reverse complement" << "\t" << endl;

                                                                    left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                    left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                    kk = kk - consective_leftSC;
                                                                    right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                    right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                    mm = mm - consective_rightSC;
                                                                    goto out_right_next;
                                                                }
                                                            } else { // No hit so check any SC pattern initiate in emitted region
                                                                // left_found = lsc, right_found = rsc
                                                                float lsc_occurance = 0;
                                                                // it = find(left_SC_base_pos.begin(), left_SC_base_pos.end(), start_lsc);
                                                                bool left_found_in_left = false;

                                                                for (it = left_SC_base_pos.begin(); it != left_SC_base_pos.end(); it++) {

                                                                    if (abs(*it - start_lsc) <= left_right_dist) {
                                                                        left_found_in_left = true;
                                                                        break;
                                                                    }
                                                                }
                                                                if (left_found_in_left) {//if (it != left_SC_base_pos.end()) {
                                                                    int pos = std::distance(left_SC_base_pos.begin(), it);
                                                                    vector <string> seq_lsc_occ;

                                                                    for (int l_found = pos; left_SC_base_pos[l_found] != 0; l_found++) {
                                                                        seq_lsc_occ.push_back(left_SC_seq[l_found]);
                                                                        lsc_occurance++;
                                                                    }

                                                                    if (lsc_occurance >= consective_SC_TH && start_lsc != previous_last_leftSC) { // lsc pattern found at the start of emitted sequence
                                                                        int pos_check;

                                                                        if (cov.at(start_lsc) > 1) // to avoid zero cov at sc position -2 for RSC none for LSC
                                                                            pos_check = start_lsc;
                                                                        else
                                                                            pos_check = start_lsc - 1;

                                                                        if (float(lsc_occurance / cov.at(pos_check - 1)) * 100 >= SC_support_TH2) { // -1 because coverage starts from 0 index not one

                                                                            string cons_lsc_occ = utils.consensus_seq_CAP3(seq_lsc_occ, toleft, PATH); //"left");

                                                                            if (cons_lsc_occ.size() >= cons_seq_size_TH) {

                                                                                string cons_left_blast = perform_blast(cons_lsc_occ, contig_seq); //boost::regex_search(contig_seq, present, regexPattern);
                                                                                vector<string> temp3;// = utils.split(cons_left_blast, "\t");
                                                                                utils.str_split(cons_left_blast, temp3, delimiter);
                                                                                int cons_left_found = atoi(temp3[0].c_str());
                                                                                int cons_left_end = atoi(temp3[2].c_str()) - cons_left_found;
                                                                                string cons_left_strand = temp3[1];
                                                                                left_st_ind = pos;
                                                                                left_en_ind = pos + lsc_occurance;
                                                                                
                                                                                if (cons_left_found != -1) {
                                                                                    
                                                                                    if (cons_left_strand == "+") {
                                                                                        cons_left_found = cons_left_found + cons_left_end; // end pos of match to insert
                                                                                        problem_identified = true; // *translocation* with left and right fragments inversion and translocation
                                                                                        sc_patterns.push_back(start_lsc);
                                                                                        sc_patterns.push_back(end_lsc);
                                                                                        sc_patterns.push_back(found_left);
                                                                                        out_file << contig_name << "\t" << "inversion_translocation_translocation" << "\t" << end_lsc << "-" << contig_seq.size() << "\t" << " at:" << found_left << " and translocate: " << start_lsc << "-" << end_lsc << "\t" << "forward" << " at :" << cons_left_found << "\t" << endl;
                                                                                        left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                                        left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                                        kk = kk - consective_leftSC;
                                                                                        right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                                        right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                                        mm = mm - consective_rightSC;
                                                                                        goto out_right_next;
                                                                                        
                                                                                    } else {// if (cons_left_strand == "-") {
                                                                                        cons_lsc_occ = utils.Rcomplement(cons_lsc_occ); // take reverse complement of LSC seq 

                                                                                        problem_identified = true; // *inversion and translocation* with left and right fragments inversion and translocation
                                                                                        sc_patterns.push_back(start_lsc);
                                                                                        sc_patterns.push_back(end_lsc);
                                                                                        sc_patterns.push_back(found_left);
                                                                                        sc_patterns.push_back(cons_left_found);
                                                                                        out_file << contig_name << "\t" << "inversion_translocation_inversion_translocation" << "\t" << end_lsc << "-" << contig_seq.size() << "\t" << "invert and translocate at:" << found_left << "\t" << start_lsc << "-" << end_lsc << "\t" << "reverse complement" << "\t" << "at :" << cons_left_found << endl;
                                                                                        left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                                        left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                                        kk = kk - consective_leftSC;
                                                                                        right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                                        right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                                        mm = mm - consective_rightSC;
                                                                                        goto out_right_next;
                                                                                    }
                                                                                }
                                                                            } // consenses check
                                                                        } //cov check
                                                                    }
                                                                    lsc_occurance = 0;
                                                                } // no RSC check because its already been checked to be not found
                                                            }
                                                            // }

                                                            if (!problem_identified) {
                                                                int st =  found_left_end; // found_left + 17jan
                                                                query_frag = contig_seq.substr(st, last_leftSC - st);

                                                                problem_identified = true; // simple inversion and translocation of query fragment, insert at found_left
                                                                sc_patterns.push_back(start_lsc);
                                                                sc_patterns.push_back(end_lsc);
                                                                sc_patterns.push_back(found_left);
                                                                out_file << contig_name << "\t" << "inversion_translocation" << "\t" << start_lsc << "-" << end_lsc << "\t" << "Reverse complement" << "\t" << "insert at:" << found_left << endl;
                                                              //  out_file << contig_name << "\t" << "inversion_translocation" << "\t" << end_lsc << "-" << contig_seq.size() << "\t" << "Reverse complement" << "\t" << "invert  and translocate at:" << found_left << endl;

                                                                left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                kk = kk - consective_leftSC;
                                                                right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                mm = mm - consective_rightSC;
                                                                goto out_right_next;
                                                            }

                                                        } // else LSC seq not found  
                                                    } // cons check
                                                    //LSC seq either not found so it can be missingSequence #### 001
                                                    if (!problem_identified) { // LSC couldn't help means it is crossed pattern with missingSequence
                                                        
                                                        for (vector<int>::iterator c = cov.begin(); c != cov.end(); c++) {
                                                            
                                                            int i = std::distance(cov.begin(), c) - 1; // 1 because cov index from 0;
                                                            
                                                            if ((i == last_leftSC && float(consective_leftSC / *c)*100 >= SC_support_TH2) || (i == last_rightSC && float(consective_rightSC / (*c) * 100 >= SC_support_TH2))) {
                                                                problem_identified = true; // missingSequence between last_rightSC and last_leftSC    
                                                                sc_patterns.push_back(last_rightSC);
                                                                sc_patterns.push_back(last_leftSC);
                                                                //*********
                                                                // find overlap between LSC and RSC and make proper sequence to put here
                                                                // RSC + Overlap(if any) + LSC
                                                                //*********

                                                                left_cons_seq = utils.consensus_seq_CAP3(left_seq, toleft, PATH); //"left");
                                                                right_cons_seq = utils.consensus_seq_CAP3(right_seq, !toleft, PATH); //"right");

                                                                if (!left_cons_seq.empty() && !right_cons_seq.empty()) {

                                                                    /**16Dec21
                                                                     filter missing sequence position using direction of mapped reads at that position*/

                                                                    bool isMissingSequence = filter_missingSeq(contig_name, last_rightSC, bam_file);

                                                                    if (isMissingSequence) {

                                                                        string RSC_LSC = overlapSC(right_cons_seq, left_cons_seq);

                                                                        if (!RSC_LSC.empty()) {
                                                                            // insert RSC_LSC at previous_last_rightSC or previous_last_leftSC
                                                                            out_file << contig_name << "\t" << "missingSequence (w overlap)" << "\t" << last_rightSC << "\t" << "forward" << "\t" << RSC_LSC << endl;
                                                                            left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                            left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                            right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                            right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                            mm = mm - consective_rightSC;
                                                                        } else {
                                                                            //insert RSC + LSC at previous_last_rightSC or previous_last_leftSC
                                                                            out_file << contig_name << "\t" << "missingSequence (w/o overlap)" << "\t" << last_rightSC << "\t" << "forward" << "\t" << right_cons_seq << left_cons_seq << endl;
                                                                            left_SC_base_pos.erase(left_SC_base_pos.begin() + left_st_ind, left_SC_base_pos.begin() + left_en_ind);
                                                                            left_SC_seq.erase(left_SC_seq.begin() + left_st_ind, left_SC_seq.begin() + left_en_ind);
                                                                            right_SC_base_pos.erase(right_SC_base_pos.begin() + (mm - consective_rightSC), right_SC_base_pos.begin() + mm);
                                                                            right_SC_seq.erase(right_SC_seq.begin() + (mm - consective_rightSC), right_SC_seq.begin() + mm);
                                                                            mm = mm - consective_rightSC;
                                                                        }
                                                                        crossed_found = true;
                                                                        goto out_right_next;
                                                                    } // filter missing sequence positions 
                                                                } // consensus check
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            left_seq.clear();
                                            last_leftSC = 0;
                                            consective_leftSC = 0;
                                        }
                                    }
                                }
                                left_seq.clear();
                                last_leftSC = 0;
                                consective_leftSC = 0;
                                // add thing of smaller and larger fragment here
                                if (!(crossed_found)) { // should check coverage here , can be self or multigeneChimera 
                                    int start, end;
                                    first_frag = contig_seq.substr(0, last_rightSC);
                                    if (last_rightSC + 1 < contig_seq.size())
                                        sec_frag = contig_seq.substr(last_rightSC + 1, contig_seq.size() - 1); // fragment after
                                    else
                                        sec_frag = "";
                                    if (first_frag.size() < sec_frag.size()) { // take smallest fragment of contig to check for chimera or local mis-assembly
                                        query_frag = first_frag;
                                        first_frag = "";
                                        start = 0;
                                        end = last_rightSC;
                                    } else { // sec fragment is smaller
                                        query_frag = sec_frag;
                                        sec_frag = "";
                                        start = last_rightSC;
                                        end = contig_seq.length() - 1;
                                    }
                                    string target = first_frag + string(query_frag.size(), 'N') + sec_frag;
                                    string blast_hit = perform_blast(query_frag, first_frag, target);
                                    vector<string> temp2;// = utils.split(blast_hit, "\t");
                                    utils.str_split(blast_hit, temp2, delimiter);
                                    blast_found = temp2[0];
                                    
                                    int target_start = atoi(temp2[1].c_str());
                                    int target_end = atoi(temp2[2].c_str());
                                    bool repeat = false;
                                    
                                    if (!sc_patterns.empty()) { // check miss-assembled regions
                                        
                                        for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                            
                                            if (abs(target_start - *sc) < repeat_TH) {
                                                
                                                for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                    
                                                    if ((abs(target_end - *sc) < repeat_TH)) {
                                                        repeat = true;
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    if (repeat) {

                                        goto out_right_next;
                                    }
                                    if (blast_found != "*") {
                                        sc_patterns.push_back(target_start);
                                        sc_patterns.push_back(target_end);
                                        problem_identified = true; // selfChimera  start and end pos
                                        sc_patterns.push_back(start);
                                        sc_patterns.push_back(end);
                                        sc_patterns.push_back(target_start);
                                        sc_patterns.push_back(target_end);
                                        
                                        if (blast_found == "+") {
                                            out_file << contig_name << "\t" << "selfChimera" << "\t" << start << "-" << end << "\t" << "Forward" << "\t" << "Discard" << endl;
                                        } else {
                                            query_frag = utils.Rcomplement(query_frag);
                                            out_file << contig_name << "\t" << "selfChimera" << "\t" << start << "-" << end << "\t" << "Reverse complement" << "\t" << "Discard" << endl;
                                        }
                                        goto out_right_next;

                                    } else {
                                        if (abs(end - start) >= ignore_short_seq) { // no cov check because its either one side of contig, not in the middle 

                                            //28oct21 avoid chimeras at the start or end positions of overlapped merged positions found in RI and merged by SCs
                                            float SC_TH_mc;

                                            if (find(left_SC_base_pos.begin(), left_SC_base_pos.end(), start - 1) != left_SC_base_pos.end()) { // if crisscross pattern
                                                SC_TH_mc = SC_support_TH2;
                                            } else // single sided pattern
                                                SC_TH_mc = SC_support_TH;
                                            // cout << cov.at(pos_check - 1) << endl;
                                            if (float(consective_rightSC / cov.at(pos_check - 1)) * 100 >= SC_TH_mc) { // 10feb22

                                            if (merged_contigs_itr == merged_contigs.end()) { 
                                                problem_identified = true;
                                                sc_patterns.push_back(start);
                                                sc_patterns.push_back(end);
                                                out_file << contig_name << "\t" << "multigeneChimera" << "\t" << start << "-" << end << "\t" << "Forward" << "\t" << "New contig" << endl;
                                                
                                            } else { //identified chimera positions must be at one side of merged points in the merged contigs

                                                if ((start < start_overlap_merged && end < start_overlap_merged) || (start > end_overlap_merged && end > end_overlap_merged)) {
                                                    problem_identified = true;
                                                    sc_patterns.push_back(start);
                                                    sc_patterns.push_back(end);
                                                    out_file << contig_name << "\t" << "multigeneChimera" << "\t" << start << "-" << end << "\t" << "Forward" << "\t" << "New contig" << endl;
                                                } // else go out
                                            }
                                            }
                                            goto out_right_next;
                                        } else {
                                            out_file << contig_name << "\t" << "miss-assembly" << "\t" << start << "-" << end << "\t" << "Forward" << "\t" << "discard" << endl;
                                            goto out_right_next;
                                        }
                                    }
                                    //}
                                }
                                //     }
                            }
                        } // consenses sequence is not empty
                        //if RSC meet certain coverage criteria
                    }

                }
out_right_next:
                ;
                problem_identified = false;
                right_seq.clear();
                last_rightSC = 0;
                consective_rightSC = 0;
                left_seq.clear();
                last_leftSC = 0;
                consective_leftSC = 0;
            }
        }
        problem_identified = false;
        right_seq.clear();
        last_rightSC = 0;
        consective_rightSC = 0;
        //  }
    }
    // LSC which are left after treating with RSC above
    if (!(left_SC_base_pos.empty())) { // no RSC only left, its single sided sc pattern, 
        int start_lsc, end_lsc;
        
        for (int mm = 0; mm < left_SC_base_pos.size(); mm++) { //for 1 left SC cluster check all right soft clip clusters ...between two clusters 0 is present to differentiate them
            
            if (left_SC_base_pos[mm] != 0) {
                last_leftSC = left_SC_base_pos[mm];
                consective_leftSC += 1;
                left_seq.push_back(left_SC_seq[mm]);
            } else {
                
                if (consective_leftSC >= consective_SC_TH && abs(last_leftSC - previous_last_leftSC) > consecutive_missAssembled_pos_dist) {
                    int pos_check;
                    
                    if (cov.at(last_leftSC - 1 ) > 1) // to avoid zero cov at sc position -2 for RSC none for LSC
                        pos_check = last_leftSC;
                    else 
                        pos_check = last_leftSC - 1;

                    if (float(consective_leftSC / cov.at(pos_check - 1)) * 100 >= SC_support_TH2) { // -1 because coverage starts from 0 index not one 
                    //    cout << "pos_checked" << endl;
                        left_cons_seq = utils.consensus_seq_CAP3(left_seq, toleft, PATH);//"left");
                        
                        if (left_cons_seq.size() >= cons_seq_size_TH) {
                            previous_last_leftSC = last_leftSC;

                            string found_left_blast = perform_blast(left_cons_seq, contig_seq); //boost::regex_search(contig_seq, present, regexPattern);
                            vector<string> temp2 ;//= utils.split(found_left_blast, "\t");
                            utils.str_split(found_left_blast, temp2, delimiter);
                            found_left = atoi(temp2[0].c_str());
                            found_left_end = atoi(temp2[2].c_str()); // 17jan - found_left; 
                            string found_left_strand = temp2[1];
                            
                            if (found_left != -1 && abs(found_left_end - last_leftSC) > left_right_dist){ //{ if (found_left != -1) { // LSC seq found in contig -> merge two points 
                                
                                if (found_left_strand == "+") {
                                    //found_left = present.position(a);
                                    //extract middle part
                                    // do blast -> self both for and RC
                                    // else find pos < lsc pos -> check any LSC at start pos -> forword?  translocation-> RC? inversion and translocation(rsc ALREADY CEHCKED)
                                    // if pos > lsc simple translocationof fragment between pos -lsc insert at end of match found  
                                    if (last_leftSC > found_left + 10) { // 16 june 2021
                                        start_lsc =  found_left_end; //found_left + // 17jan
                                        end_lsc = last_leftSC - 1;
                                        query_frag = contig_seq.substr(start_lsc, end_lsc - start_lsc); // middle suspicious fragment
                                        first_frag = contig_seq.substr(0, start_lsc - 1); // fragment before suspicious fragment
                                        
                                        if (end_lsc + 1 < contig_seq.size())
                                            sec_frag = contig_seq.substr(end_lsc + 1, contig_seq.size() - 1); // fragment after
                                        else
                                            sec_frag = "";
                                    } else { // if (last_leftSC <= found_left) { // SC seq fond on left side of contig{
                                          // 28/12/19 edited
//                                        end_lsc = found_left, start_lsc = last_leftSC;
//                                        query_frag = contig_seq.substr(start_lsc, end_lsc - start_lsc); // middle suspicious fragment
//                                        problem_identified = true; // translocationcut from start to end and put at found_left + cons_seq.size();
//                                        sc_patterns.push_back(start_lsc);
//                                        sc_patterns.push_back(end_lsc);
//                                        sc_patterns.push_back(found_left);
//                                        out_file << contig_name << "\t" << "translocation" << "\t" << start_lsc << "-" << end_lsc << "\t" << "forward" << "\t" << "translocate at:" << found_left + left_cons_seq.size() - 1 << endl;
                                        goto out_left_next;
                                    }

                                    //come here only when last_leftSC  > found_left
                                    //start******************************************************************
                                    string target = first_frag + string(query_frag.size(), 'N') + sec_frag;
                                    string blast_hit = perform_blast(query_frag, target, "");
                                    vector<string> temp2;// = utils.split(blast_hit, "\t");
                                    utils.str_split(blast_hit, temp2, delimiter);
                                    blast_found = temp2[0];
                                    
                                    int target_start = atoi(temp2[1].c_str());
                                    int target_end = atoi(temp2[2].c_str());

                                    bool repeat = false;
                                    if (!sc_patterns.empty()) { // check miss-assembled regions
                                        
                                        for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                            
                                            if (abs(target_start - *sc) < repeat_TH) {
                                                
                                                for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                    
                                                    if ((abs(target_end - *sc) < repeat_TH)) {
                                                        repeat = true;
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    if (repeat) {
                                        //                                            sc_patterns.push_back(target_start);
                                        //                                            sc_patterns.push_back(target_end);
                                        goto out_left_next;
                                        
                                    } else if (blast_found != "*") {
                                        sc_patterns.push_back(target_start);
                                        sc_patterns.push_back(target_end);
                                        
                                        if (blast_found == "+") {
                                            problem_identified = true; //selfChimera -> merge first and sec fragment
                                            sc_patterns.push_back(start_lsc);
                                            sc_patterns.push_back(end_lsc);
                                            out_file << contig_name << "\t" << "selfChimera" << "\t" << start_lsc << "-" << end_lsc << "\t" << "forward" << "\t" << "Discard" << endl;
                                            goto out_left_next;

                                        } else { // RC of query seq
                                            query_frag = utils.Rcomplement(query_frag);
                                            // blast_found = perform_blast(query_frag, first_frag, sec_frag);
                                            problem_identified = true; //selfChimera -> merge first and sec fragment
                                            sc_patterns.push_back(start_lsc);
                                            sc_patterns.push_back(end_lsc);
                                            out_file << contig_name << "\t" << "selfChimera" << "\t" << start_lsc << "-" << end_lsc << "\t" << "Reverse complement" << "\t" << "Discard" << endl;
                                            goto out_left_next;
                                        }
                                    } else { // check any SC pattern initiate from emitted region (query fragment)
                                        // left_found = lsc, right_found = rsc
                                        float lsc_occurance = 0;
                                       // it = find(left_SC_base_pos.begin(), left_SC_base_pos.end(), found_left);

                                        bool left_found_in_left = false;

                                        for (it = left_SC_base_pos.begin(); it != left_SC_base_pos.end(); it++) {

                                            if (abs(*it - found_left) <= left_right_dist) {
                                                left_found_in_left = true;
                                                break;
                                            }
                                        }

                                        if (left_found_in_left) { //(it != left_SC_base_pos.end()) {
                                            int pos = std::distance(left_SC_base_pos.begin(), it);
                                            vector <string> seq_lsc_occ;

                                            for (int l_found = pos; left_SC_base_pos[l_found] != 0; l_found++) {
                                                seq_lsc_occ.push_back(left_SC_seq[l_found]);
                                                lsc_occurance++;
                                            }
                                            if (lsc_occurance >= consective_SC_TH && start_lsc != previous_last_leftSC) { // lsc pattern found at the start of emitted sequence
                                                
                                                int pos_check;
                                                
                                                if (cov.at(start_lsc) > 1) // to avoid zero cov at sc position -2 for RSC none for LSC
                                                    pos_check = start_lsc ;
                                                else 
                                                    pos_check = start_lsc - 1;
                                                
                                                if (float(lsc_occurance / cov.at(pos_check - 1)) * 100 >= SC_support_TH2) { // -1 because coverage starts from 0 index not one
                                                    string cons_lsc_occ = utils.consensus_seq_CAP3(seq_lsc_occ, toleft, PATH);//"left");
                                                    
                                                    if (cons_lsc_occ.size() >= cons_seq_size_TH) {

                                                        string cons_left_blast = perform_blast(cons_lsc_occ, contig_seq); //boost::regex_search(contig_seq, present, regexPattern);
                                                        vector<string> temp3;// = utils.split(cons_left_blast, "\t");
                                                        utils.str_split(cons_left_blast, temp3, delimiter);
                                                        int cons_left_found = atoi(temp3[0].c_str());
                                                        int cons_left_end = atoi(temp3[2].c_str()) - cons_left_found;
                                                        string con_left_strand = temp3[1];
                                                        
                                                        if (cons_left_found != -1) {
                                                            
                                                            if (con_left_strand == "+") {
                                                                //cons_left_found = cons_left_found + cons_left_end; // end pos of match to insert
                                                                problem_identified = true; // *translocation* insert query after cons_left_found and merge first/sec frags
                                                                sc_patterns.push_back(start_lsc);
                                                                sc_patterns.push_back(end_lsc);
                                                                sc_patterns.push_back(cons_left_found);
                                                                out_file << contig_name << "\t" << "translocation" << "\t" << start_lsc << "-" << end_lsc << "\t" << "forward" << "\t" << "translocate at:" << cons_left_found + cons_lsc_occ.size() - 1 << endl;
                                                                goto out_left_next;
                                                                
                                                            } else{ // if (con_left_strand == "-") {
                                                                cons_lsc_occ = utils.Rcomplement(cons_lsc_occ); // take reverse complement of LSC seq 
                                                                problem_identified = true; // *inversion and translocation* insert query after cons_left_found and merge first/sec frags
                                                                sc_patterns.push_back(start_lsc);
                                                                sc_patterns.push_back(end_lsc);
                                                                sc_patterns.push_back(cons_left_found);
                                                                out_file << contig_name << "\t" << " inversion_translocation" << "\t" << start_lsc << "-" << end_lsc << "\t" << "Reverse complement" << "\t" << "translocate at:" << cons_left_found << endl;
                                                                goto out_left_next;
                                                            }
                                                        }
                                                    }// concensus check
                                                }
                                            }
                                            lsc_occurance = 0;
                                        }
                                        // if problem not identified yet don't check right as its already been checked above
                                        //3/12/2019
                                        if (!problem_identified) { // before moving towards multigeneChimera go for translocatin
                                            first_frag = contig_seq.substr(0, last_leftSC);
                                            int start, end;
                                            if (last_leftSC + 1 < contig_seq.size())
                                            sec_frag = contig_seq.substr(last_leftSC, contig_seq.size() - 1); // fragment after
                                        else
                                            sec_frag = "";

                                            if (first_frag.size() < sec_frag.size()) { // take smallest fragment of contig to check for chimera or local mis-assembly
                                                query_frag = first_frag;
                                                first_frag = "";
                                                start = 0;
                                                end = last_leftSC;
                                                string blast_hit = perform_blast(query_frag, sec_frag, "");
                                                
                                            } else { // sec fragment is smaller
                                                query_frag = sec_frag;
                                                sec_frag = "";
                                                start = last_leftSC;
                                                end = contig_seq.length() - 1;
                                                string blast_hit = perform_blast(query_frag, first_frag, "");
                                            }
                                            
                                            vector<string> temp2;// = utils.split(blast_hit, "\t");
                                            utils.str_split(blast_hit, temp2, delimiter);
                                            blast_found = temp2[0];
                                            int target_start = atoi(temp2[1].c_str());
                                            int target_end = atoi(temp2[2].c_str());
                                            bool repeat = false;
                                            
                                            if (!sc_patterns.empty()) { // check miss-assembled regions
                                            
                                                for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                    
                                                    if (abs(target_start - *sc) < repeat_TH) {
                                                    
                                                        for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                        
                                                            if ((abs(target_end - *sc) < repeat_TH)) {
                                                                repeat = true;
                                                                break;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            if (repeat) {
                                                goto out_left_next;
                                            }
                                            
                                            if (blast_found != "*") {
                                                sc_patterns.push_back(target_start);
                                                sc_patterns.push_back(target_end);
                                                problem_identified = true; // selfChimera  start and end pos
                                                sc_patterns.push_back(start);
                                                sc_patterns.push_back(end);
                                                sc_patterns.push_back(target_start);
                                                sc_patterns.push_back(target_end);
                                                
                                                if (blast_found == "+") {
                                                    out_file << contig_name << "\t" << "selfChimera" << "\t" << start << "-" << end << "\t" << "Forward" << "\t" << "Discard" << endl;
                                                } else {
                                                    query_frag = utils.Rcomplement(query_frag);
                                                    out_file << contig_name << "\t" << "selfChimera" << "\t" << start << "-" << end << "\t" << "Reverse complement" << "\t" << "Discard" << endl;
                                                }
                                                goto out_left_next;

                                            } else {
                                                int sum = 0;
                                                for (int aa = start_lsc; aa <= end_lsc; aa++) {
                                                    sum = sum + cov[aa];
                                                }
                                                int query_frag_avg = float(sum) / float(query_frag.size());

                                                if (abs(end_lsc - start_lsc) >= ignore_short_seq && query_frag_avg >= query_frag_avg_TH) {
                                                    
                                                    if (!sc_patterns.empty()) { // check miss-assembled regions
                                                        
                                                        for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                           
                                                            if (abs(start_lsc - *sc) < repeat_TH) {
                                                                
                                                                for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                   
                                                                    if ((abs(end_lsc - *sc) < repeat_TH)) {
                                                                        repeat = true;
                                                                        break;
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                            if (repeat) {
                                                            goto out_left_next;
                                                    }

                                                    //28oct21 avoid chimeras at the start or end positions of overlapped merged positions found in RI and merged by SCs
                                                    float SC_TH_mc;
                                                        
                                                    if (find(right_SC_base_pos.begin(), right_SC_base_pos.end(), last_leftSC  + 1) != right_SC_base_pos.end()) { // if crisscross pattern
                                                        SC_TH_mc = SC_support_TH2;
                                                    } else // single sided pattern
                                                        SC_TH_mc = SC_support_TH;
                                                        //cout <<float(consective_leftSC / cov.at(pos_check - 1)) << endl;
                                                    if (float(consective_leftSC / cov.at(pos_check - 1)) * 100 >= SC_TH_mc) { // 10feb22

                                                    if (merged_contigs_itr == merged_contigs.end()) {

                                                        problem_identified = true; //trasns chimera 
                                                        sc_patterns.push_back(start_lsc);
                                                        sc_patterns.push_back(end_lsc);
                                                        out_file << contig_name << "\t" << "multigeneChimera" << "\t" << start_lsc << "-" << end_lsc << "\t" << "forward" << "\t" << "New contig" << endl;

                                                    } else { //identified chimera positions must be at one side of merged points in the merged contigs
                                                        if ((start_lsc < start_overlap_merged && end_lsc < start_overlap_merged) || (start_lsc > end_overlap_merged && end_lsc > end_overlap_merged)) {
                                                            problem_identified = true; //trasns chimera 
                                                            sc_patterns.push_back(start_lsc);
                                                            sc_patterns.push_back(end_lsc);
                                                            out_file << contig_name << "\t" << "multigeneChimera" << "\t" << start_lsc << "-" << end_lsc << "\t" << "forward" << "\t" << "New contig" << endl;
                                                        } // else go out
                                                    }
                                                }
                                                    goto out_left_next;

                                                } else if (query_frag.size() < ignore_short_seq) {
                                                    problem_identified = true; // merge start and end discard middle part
                                                    sc_patterns.push_back(start_lsc);
                                                    sc_patterns.push_back(end_lsc);
                                                    out_file << contig_name << "\t" << "unsupported_insertion" << "\t" << start_lsc << "-" << end_lsc << "\t" << "forward" << "\t" << "Discard" << endl;
                                                    goto out_left_next;

                                                } else {
                                                        sc_patterns.push_back(last_leftSC);
                                                        sc_patterns.push_back(contig_seq.size() - 1);
                                                        sc_patterns.push_back(found_left + left_cons_seq.size() - 1);
                                                        query_frag = contig_seq.substr(start_lsc, end_lsc - start_lsc); // join two points and emit middle part // when last_leftSC > found_left
                                                        out_file << contig_name << "\t" << "translocation" << "\t" << last_leftSC << "-" << contig_seq.size() - 1 << "\t" << "forward" << "\t" << "translocate at:" << found_left + left_cons_seq.size() - 1 << endl;
                                                        problem_identified = true;
                                                goto out_left_next;
                                                //apply length and coverage threshold
                                        }    }
                                        
                                        }
                                    }
                                    //  }
                                    //end ******************************************************************************

//                                    if (!problem_identified) { // check if query fragment meet criteria of multigeneChimera because LSC seq found somewhere in contig while checking for crossed pattern
//                                        query_frag = contig_seq.substr(start_lsc, end_lsc - start_lsc); // join two points and emit middle part // when last_leftSC > found_left
//
//                                        if (query_frag.size() >= 100) {
//                                            if (!sc_patterns.empty()) { // check miss-assembled regions
//                                                for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
//                                                    if (abs(start_lsc - *sc) < 100) {
//                                                        for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
//                                                            if ((abs(end_lsc - *sc) < 100)) {
//                                                                repeat = true;
//                                                                break;
//                                                            }
//                                                        }
//                                                    }
//                                                }
//                                            }
//                                            if (repeat == true) {
//                                                goto out_left_next;
//                                            }
//                                            problem_identified = true; //trasns chimera 
//                                            sc_patterns.push_back(start_lsc);
//                                            sc_patterns.push_back(end_lsc);
//                                            out_file << contig_name << "\t" << "multigeneChimera" << "\t" << start_lsc << "-" << end_lsc << "\t" << "forward" << "\t" << "New contig" << endl;
//                                            goto out_left_next;
//                                        } else if (query_frag.size() < ignore_short_seq){  // 17 may 2021 bug fixed {
//                                            problem_identified = true; // merge start and end discard middle part
//                                            sc_patterns.push_back(start_lsc);
//                                            sc_patterns.push_back(end_lsc);
//                                            out_file << contig_name << "\t" << "unsupported_insertion" << "\t" << start_lsc << "-" << end_lsc << "\t" << "forward" << "\t" << "Discard" << endl;
//                                            goto out_left_next;
//                                        }
//                                    }

                                } else if (found_left_strand == "-") { // for LSC seq not found in contig search RC of LSC
                                    left_cons_seq = utils.Rcomplement(left_cons_seq);

                                    //found_left = present.position(a);
                                    if (last_leftSC > found_left + 10) { // 16 june 2021
                                        start_lsc = found_left + left_cons_seq.size() + 1;
                                        end_lsc = last_leftSC;
                                        query_frag = contig_seq.substr(start_lsc, end_lsc - start_lsc); // middle suspicious fragment
                                        first_frag = contig_seq.substr(0, start_lsc - 1); // fragment before suspicious fragment          
                                        
                                        if (end_lsc + 1 < contig_seq.size())
                                            sec_frag = contig_seq.substr(end_lsc + 1, contig_seq.size() - 1); // fragment after
                                        else
                                            sec_frag = "";
                                    } else{// if (last_leftSC < found_left) { // SC seq fond on right side of contig{
                                        sc_patterns.push_back(start_lsc);
                                        sc_patterns.push_back(end_lsc);
                                        end_lsc = found_left, start_lsc = last_leftSC;
                                        query_frag = contig_seq.substr(start_lsc, end_lsc - start_lsc); // middle suspicious fragment
                                        problem_identified = true; // translocationcut from start to end and invert at its own place
                                        
                                        out_file << contig_name << "\t" << "inversion" << "\t" << start_lsc << "-" << end_lsc << "\t" << "Reverse complement" << "\t" << "invert at its own place" << endl;
                                        goto out_left_next;
                                    }
                                    string target = first_frag + string(query_frag.size(), 'N') + sec_frag;
                                    string blast_hit = perform_blast(query_frag, target, "");
                                    
                                    vector<string> temp2;// = utils.split(blast_hit, "\t");
                                    
                                    utils.str_split(blast_hit, temp2, delimiter);
                                    blast_found = temp2[0];
                                    
                                    int target_start = atoi(temp2[1].c_str());
                                    int target_end = atoi(temp2[2].c_str());
                                    bool repeat = false;
                                    
                                    if (!sc_patterns.empty()) { // check miss-assembled regions
                                        
                                        for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                            
                                            if (abs(target_start - *sc) < repeat_TH) {
                                                
                                                for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                    
                                                    if ((abs(target_end - *sc) < repeat_TH)) {
                                                        repeat = true;
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    if (repeat) {
                                        //                                            sc_patterns.push_back(target_start);
                                        //                                            sc_patterns.push_back(target_end);
                                        goto out_left_next;
                                        
                                    } else if (blast_found != "*") {
                                        sc_patterns.push_back(target_start);
                                        sc_patterns.push_back(target_end);
                                        problem_identified = true; //selfChimera with inversion and translocation of left frag at start of match
                                        sc_patterns.push_back(start_lsc);
                                        sc_patterns.push_back(end_lsc);
                                        
                                        if (blast_found == "+") {
                                            out_file << contig_name << "\t" << "Self chimera" << "\t" << start_lsc << "-" << end_lsc << "\t" << "forward" << "\t" << "Discard" << endl;

                                        } else { // RC of query seq
                                            query_frag = utils.Rcomplement(query_frag);
                                            out_file << contig_name << "\t" << "Self chimera" << "\t" << start_lsc << "-" << end_lsc << "\t" << "Reverse complement" << "\t" << "Discard" << endl;
                                        }
                                        goto out_left_next;
                                        
                                    } else { // check any SC pattern initiate in emitted region
                                        // left_found = lsc, right_found = rsc
                                        if(query_frag.size() >= ignore_short_seq){

                                            float lsc_occurance = 0;
                                            //it = find(left_SC_base_pos.begin(), left_SC_base_pos.end(), found_left);

                                            bool left_found_in_left = false;

                                            for (it = left_SC_base_pos.begin(); it != left_SC_base_pos.end(); it++) {

                                                if (abs(*it - found_left) <= left_right_dist) {
                                                    left_found_in_left = true;
                                                    break;
                                                }
                                            }

                                            if(left_found_in_left){ //if (it != left_SC_base_pos.end()) {
                                                int pos = std::distance(left_SC_base_pos.begin(), it);
                                                vector <string> seq_lsc_occ;

                                                for (int l_found = pos; left_SC_base_pos[l_found] != 0; l_found++) {
                                                    seq_lsc_occ.push_back(left_SC_seq[l_found]);
                                                    lsc_occurance++;
                                                }

                                            if (lsc_occurance >= consective_SC_TH && start_lsc != previous_last_leftSC) { // lsc pattern found at the start of emitted sequence
                                                int pos_check;
                                                if (cov.at(start_lsc ) > 1) // to avoid zero cov at sc position -2 for RSC none for LSC
                                                    pos_check = start_lsc ;
                                                else pos_check = start_lsc - 1;
                                                
                                                if (float(lsc_occurance / cov.at(pos_check - 1)) * 100 >= SC_support_TH2) { // -1 because coverage starts from 0 index not one
                                                    string cons_lsc_occ = utils.consensus_seq_CAP3(seq_lsc_occ, toleft, PATH);//"left");
                                                    
                                                    if (cons_lsc_occ.size() >= cons_seq_size_TH) {

                                                        string cons_left_blast = perform_blast(cons_lsc_occ, contig_seq);
                                                        vector<string> temp3 = utils.split(cons_left_blast, "\t");
                                                        
                                                        utils.str_split(cons_left_blast, temp3, delimiter);
                                                        
                                                        int cons_left_found = atoi(temp3[0].c_str());
                                                        int cons_left_end = atoi(temp3[2].c_str()) - cons_left_found;
                                                        string con_left_strand = temp3[1];
                                                        
                                                        if (cons_left_found != -1) {
                                                            
                                                            if (con_left_strand == "+") {
                                                                cons_left_found = cons_left_found + cons_left_end; // end pos of match to insert
                                                                problem_identified = true; // *translocation* with left and right fragments inversion and translocation
                                                                sc_patterns.push_back(start_lsc);
                                                                sc_patterns.push_back(end_lsc);
                                                                sc_patterns.push_back(found_left);
                                                                sc_patterns.push_back(cons_left_found);
                                                                out_file << contig_name << "\t" << "inversion_translocation_translocation" << "\t" << end_lsc << "-" << contig_seq.size() << "\t" << "invert and translocate at:" << found_left << "\t" << start_lsc << "-" << end_lsc << "\t" << "forward" << " at :" << cons_left_found << endl;

                                                                goto out_left_next;
                                                                
                                                            } else{ // if (con_left_strand == "-") {
                                                                cons_lsc_occ = utils.Rcomplement(cons_lsc_occ); // take reverse complement of LSC seq 
                                                                problem_identified = true; // *inversion and translocation* with left and right fragments inversion and translocation
                                                                sc_patterns.push_back(start_lsc);
                                                                sc_patterns.push_back(end_lsc);
                                                                sc_patterns.push_back(found_left);
                                                                sc_patterns.push_back(cons_left_found);
                                                                out_file << contig_name << "\t" << "inversion_translocation_inversion_translocation" << "\t" << end_lsc << "-" << contig_seq.size() << "\t" << "invert and translocate at:" << found_left << "\t" << start_lsc << "-" << end_lsc << "\t" << "reverse complement" << "\t" << " at :" << cons_left_found << endl;
                                                                goto out_left_next;
                                                            }
                                                        }
                                                    } // consesnsus check
                                                } // cov check
                                            }
                                            lsc_occurance = 0;
                                        } // no RSC check because its already been checked to be not found
                                    }
                                        else{
                                        out_file << contig_name << "\t" << "unsupported_insertion" << "\t" << start_lsc << "-" << end_lsc << "\t" << "forward" << "\t" << "Discard" << endl;
                                        }
                                    }
                                    // }

                                    if (!problem_identified) {
                                        
                                        int st = found_left + left_cons_seq.size();
                                        query_frag = contig_seq.substr(st, last_leftSC - st);
                                        problem_identified = true; // simple inversion and translocationof query fragment, insert at found_left 
                                        sc_patterns.push_back(start_lsc);
                                        sc_patterns.push_back(end_lsc);
                                        sc_patterns.push_back(found_left);
                                        //8/1/202
                                        out_file << contig_name << "\t" << "inversion_translocation" << "\t" << start_lsc << "-" << end_lsc << "\t" << "Reverse complement" << "\t" << "insert at:" << found_left << endl;
                                        //out_file << contig_name << "\t" << "inversion_translocation" << "\t" << end_lsc << "-" << contig_seq.size() << "\t" << "Reverse complement" << "\t" << "invert  and translocate at:" << found_left << endl;
                                        goto out_left_next;
                                    }

                                } // 
                            }
                            if (!problem_identified) { // check for cross pattern

                                for (int mm = 0; mm < right_SC_base_pos.size(); mm++) {
                                    //   cout<<right_SC_base_pos[kk]<<endl;
                                    if (right_SC_base_pos[mm] != 0) {
                                        last_rightSC = right_SC_base_pos[mm];
                                        consective_rightSC += 1;
                                        right_seq.push_back(right_SC_seq[mm]);

                                    } else {
                                        if (consective_rightSC >= consective_SC_TH)// && max_consec_rightSC < consective_rightSC) //max cluster supported by  >4 SC reads
                                        {
                                            right_cons_seq = utils.consensus_seq_CAP3(right_seq, !toleft, PATH);//"right");

                                            if (right_cons_seq.size() >= cons_seq_size_TH) {
                                                // to cover missing seq like TRINITY_DN20900_c0_g2
                                                if ((abs(last_leftSC - last_rightSC) <= left_right_dist) || (last_leftSC - last_rightSC >= 0 && last_leftSC < last_rightSC + right_cons_seq.length())) {

                                                    string found_right_blast = perform_blast(right_cons_seq, contig_seq);
                                                    vector<string> temp2; // = utils.split(found_right_blast, "\t");

                                                    utils.str_split(found_right_blast, temp2, delimiter);
                                                    found_right = atoi(temp2[0].c_str());

                                                    found_right_end = atoi(temp2[2].c_str()); // 17 Jan - found_right;
                                                    string found_right_strand = temp2[1];
                                                    
                                                    if (found_right != -1) {

                                                        if (found_right_strand == "+") {
                                                            int start_rsc, end_rsc;
                                                            //found_right = present.position(a);
                                                            //extract middle part
                                                            // do blast -> self both for and RC
                                                            // else find pos < lsc pos -> check any LSC at start pos -> forword?  translocation-> RC? inversion and translocation(rsc ALREADY CEHCKED)
                                                            // if pos > lsc simple translocationof fragment between pos -lsc insert at end of match found  
                                                            if (last_rightSC < found_right - 10) {  //16 june 2021
                                                                start_rsc = last_rightSC;
                                                                end_rsc = found_right;
                                                                
                                                                query_frag = contig_seq.substr(start_rsc, end_rsc - start_rsc); // middle suspicious fragment
                                                                first_frag = contig_seq.substr(0, start_rsc - 1); // fragment before suspicious fragment
                                                                
                                                                if (end_rsc + 1 < contig_seq.size())
                                                                    sec_frag = contig_seq.substr(end_rsc + 1, contig_seq.size() - 1); // fragment after
                                                                else
                                                                    sec_frag = "";
                                                            } else { // SC seq fond on right side of contig{
//                                                                start_rsc = found_right + right_cons_seq.size() + 1, end_rsc = last_rightSC;
//                                                                query_frag = contig_seq.substr(start_rsc, end_rsc); // middle suspicious fragment
//                                                                if(query_frag.size() < ignore_short_seq)
//                                                                    out_file << contig_name << "\t" << "unsupported_insertion" << "\t" << start_lsc << "-" << end_lsc << "\t" << "forward" << "\t" << "Discard" << endl;
//                                                                else{
//                                                                    
//                                                                problem_identified = true; // translocationcut from start to end and put at found_right + cons_seq.size();
//                                                                sc_patterns.push_back(start_rsc);
//                                                                sc_patterns.push_back(end_rsc);
//                                                                sc_patterns.push_back(found_right);
//                                                                out_file << contig_name << "\t" << "translocation" << "\t" << start_rsc << "-" <<  end_rsc << "\t" << "forward" << "\t" << "translocate at:" << found_right << endl;
//                                                                goto out_left_next;
//                                                            }
                                                                goto out_left_next;
                                                            }
                                                            //come here only when last_rightSC  > found_right
                                                            //start******************************************************************

                                                            string target = first_frag + string(query_frag.size(), 'N') + sec_frag;
                                                            string blast_hit = perform_blast(query_frag, target, "");
                                                            
                                                            vector<string> temp2; // = utils.split(blast_hit, "\t");
                                                            utils.str_split(blast_hit, temp2, delimiter);
                                                            blast_found = temp2[0];
                                                            
                                                            int target_start = atoi(temp2[1].c_str());
                                                            int target_end = atoi(temp2[2].c_str());
                                                            bool repeat = false;
                                                            
                                                            if (!sc_patterns.empty()) { // check miss-assembled regions
                                                                
                                                                for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                    
                                                                    if (abs(target_start - *sc) < repeat_TH) {
                                                                        
                                                                        for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                            
                                                                            if ((abs(target_end - *sc) < repeat_TH)) {
                                                                                repeat = true;
                                                                                break;
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                            if (repeat) {
                                                                //                                                                    sc_patterns.push_back(target_start);
                                                                //                                                                    sc_patterns.push_back(target_end);
                                                                goto out_left_next;
                                                            } else if (blast_found != "*") {
                                                                sc_patterns.push_back(target_start);
                                                                sc_patterns.push_back(target_end);
                                                                
                                                                if (blast_found == "+") {
                                                                    
                                                                    problem_identified = true; //selfChimera -> merge first and sec fragment
                                                                    sc_patterns.push_back(start_rsc);
                                                                    sc_patterns.push_back(end_rsc);
                                                                    out_file << contig_name << "\t" << "selfChimera" << "\t" << last_rightSC << "-" << found_right << "\t" << "forward" << "\t" << "Discard" << endl;
                                                                    goto out_left_next;

                                                                } else { // RC of query seq
                                                                    query_frag = utils.Rcomplement(query_frag);
                                                                    // blast_found = perform_blast(query_frag, first_frag, sec_frag);
                                                                    problem_identified = true; //selfChimera -> merge first and sec fragment
                                                                    sc_patterns.push_back(start_rsc);
                                                                    sc_patterns.push_back(end_rsc);
                                                                    out_file << contig_name << "\t" << "selfChimera" << "\t" << last_rightSC << "-" << found_right << "\t" << "Reverse complement" << "\t" << "Discard" << endl;
                                                                    goto out_left_next;
                                                                }
                                                            } else { // check any SC pattern initiate from emitted region (query fragment)

                                                                int sum = 0;
//                                                                for (int aa = start_rsc - 110; aa <= start_rsc - 10; aa++) {
//                                                                    sum = sum + cov[aa];
//                                                                }
//                                                                int left_avg_cov = sum / 100;
//                                                                sum = 0;
//
//                                                                for (int aa = end_rsc + 10; aa <= end_rsc + 110; aa++) {
//                                                                    sum = sum + cov[aa];
//                                                                }
//                                                                int right_avg_cov = sum / 100;
//                                                                sum = 0;

                                                                for (int aa = start_rsc; aa <= end_rsc; aa++) {
                                                                    sum = sum + cov[aa];
                                                                }
                                                                int query_frag_avg = float(sum) / float(query_frag.size());

                                                                if (abs(end_rsc - start_rsc) >= ignore_short_seq && query_frag_avg >= query_frag_avg_TH) {
                                                                    
                                                                    if (!sc_patterns.empty()) { // check miss-assembled regions
                                                                        for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                            
                                                                            if (abs(start_rsc - *sc) < repeat_TH) {
                                                                                
                                                                                for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {

                                                                                    if ((abs(end_rsc - *sc) < repeat_TH)) {
                                                                                        repeat = true;
                                                                                        break;
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                    if (repeat) {
                                                                        goto out_left_next;
                                                                    }
                                                                    //28oct21 avoid chimeras at the start or end positions of overlapped merged positions found in RI and merged by SCs
                                                                    float SC_TH_mc;

                                                                    if (find(left_SC_base_pos.begin(), left_SC_base_pos.end(), last_rightSC  - 1) != left_SC_base_pos.end()) { // if crisscross pattern
                                                                        SC_TH_mc = SC_support_TH2;
                                                                    } else // single sided pattern
                                                                        SC_TH_mc = SC_support_TH;
                                                                            //cout <<float(consective_rightSC / cov.at(pos_check - 1)) << endl;
                                                                    if (float(consective_rightSC / cov.at(pos_check - 1)) * 100 >= SC_TH_mc) { // 10feb22

                                                                    if (merged_contigs_itr == merged_contigs.end()) {

                                                                        problem_identified = true;
                                                                        sc_patterns.push_back(start_rsc);
                                                                        sc_patterns.push_back(end_rsc);
                                                                        out_file << contig_name << "\t" << "multigeneChimera" << "\t" << start_rsc << "-" << end_rsc << "\t" << "forward" << "\t" << "New contig" << endl;
                                                                        mm = mm - consective_rightSC;

                                                                    } else {//identified chimera positions must be at one side of merged points in the merged contigs

                                                                        if ((start_rsc < start_overlap_merged && end_rsc < start_overlap_merged) || (start_rsc > end_overlap_merged && end_rsc > end_overlap_merged)) {

                                                                            problem_identified = true;
                                                                            sc_patterns.push_back(start_rsc);
                                                                            sc_patterns.push_back(end_rsc);
                                                                            out_file << contig_name << "\t" << "multigeneChimera" << "\t" << start_rsc << "-" << end_rsc << "\t" << "forward" << "\t" << "New contig" << endl;
                                                                            mm = mm - consective_rightSC;
                                                                        } // else go out
                                                                    }
                                                                    }
                                                                    goto out_left_next;
                                                                } else if (query_frag.size() < ignore_short_seq) {
                                                                    problem_identified = true;
                                                                    out_file << contig_name << "\t" << "unsupported_insertion" << "\t" << start_rsc << "-" << end_rsc << "\t" << "forward" << "\t" << "Discard" << endl;
                                                                    goto out_left_next;
                                                                } else  {
                                                                     problem_identified = true;
                                                                    out_file << contig_name << "\t" << "translocation" << "\t" << 0 << "-" << last_rightSC << "\t" << "forward" << "\t" << "translocate at:" << found_right << endl;
                                                                    goto out_left_next;

                                                                }
                                                            }
                                                        } else{ // if (found_right_strand == "-") { // for LSC seq not found in contig search RC of LSC
                                                            right_cons_seq = utils.Rcomplement(right_cons_seq);
                                                            //found_right = present.position(a);
                                                            if (last_rightSC < found_right - 10) { //16june2021
                                                              
                                                                int start_rsc = last_rightSC;
                                                                int end_rsc = found_right;
                                                                query_frag = contig_seq.substr(start_rsc, end_rsc - start_rsc); // middle suspicious fragment
                                                                
                                                                first_frag = contig_seq.substr(0, start_rsc - 1); // fragment before suspicious fragment 
                                                                
                                                                if (end_rsc + 1 < contig_seq.size())
                                                                    sec_frag = contig_seq.substr(end_rsc + 1, contig_seq.size() - 1); // fragment after
                                                                else
                                                                    sec_frag = "";
                                                                
                                                            } else{// if (last_rightSC > found_right) { // SC seq fond on right side of contig{
                                                                
                                                                int end_rsc = last_rightSC, start_rsc =  found_right_end; //found_right + 17jan
                                                                //query_frag = contig_seq.substr(start_rsc, end_rsc - start_rsc); // middle suspicious fragment
                                                                problem_identified = true; // translocationcut from start to end and invert at its own place
                                                                out_file << contig_name << "\t" << "inversion" << "\t" << start_rsc << "-" << end_rsc << "\t" << "Reverse complement" << "\t" << "invert at its own place" << endl;
                                                                goto out_left_next;
                                                            }

                                                            string target = first_frag + string(query_frag.size(), 'N') + sec_frag; // extract suspicious frag and replace with Ns to maintain index number
                                                            
                                                            string blast_hit = perform_blast(query_frag, target, "");
                                                            
                                                            vector<string> temp2;// = utils.split(blast_hit, "\t");
                                                            utils.str_split(blast_hit, temp2, delimiter);
                                                            
                                                            blast_found = temp2[0];
                                                            int target_start = atoi(temp2[1].c_str());
                                                            int target_end = atoi(temp2[2].c_str());
                                                            bool repeat = false;
                                                            
                                                            if (!sc_patterns.empty()) { // check miss-assembled regions
                                                                
                                                                for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                    
                                                                    if (abs(target_start - *sc) < repeat_TH) {
                                                                        
                                                                        for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                                           
                                                                            if ((abs(target_end - *sc) < repeat_TH)) {
                                                                                repeat = true;
                                                                                break;
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                            if (repeat) {
                               
                                                                goto out_left_next;
                                                            } else if (blast_found != "*") { // hit found
                                                                sc_patterns.push_back(target_start);
                                                                sc_patterns.push_back(target_end);
                                                                problem_identified = true; //selfChimera with inversion and translocationof right frag at start of match
                                                                sc_patterns.push_back(found_right);
                                                                sc_patterns.push_back(last_rightSC);
                                                                
                                                                if (blast_found == "+") {
                                                                    out_file << contig_name << "\t" << "Self chimera" << "\t" << last_rightSC << "-" << found_right << "\t" << "forward" << "\t" << "Discard" << endl;

                                                                } else { // RC of query seq
                                                                    query_frag = utils.Rcomplement(query_frag);
                                                                    out_file << contig_name << "\t" << "Self chimera" << "\t" << last_rightSC << "-" << found_right << "\t" << "Reverse complement" << "\t" << "Discard" << endl;
                                                                }
                                                                goto out_left_next;
                                                            } else { // hit not found, check any SC pattern initiate in emitted region
                                                               // int st = found_right + right_cons_seq.size(); //commented 17 jan 21
                                                                //query_frag = contig_seq.substr(st, last_rightSC - st);
                                                                problem_identified = true; // simple inversion and translocationof query fragment, insert at found_right 
                                                                sc_patterns.push_back(last_rightSC);
                                                                sc_patterns.push_back(found_right);
                                                                sc_patterns.push_back(last_leftSC);
                                                                out_file << contig_name << "\t" << "inversion_translocation" << "\t" << 0 << "-" << last_rightSC << "\t" << "Reverse complement" << "\t" << "invert  and translocate at:" <<  found_right_end  << endl;
                                                                goto out_left_next;

                                                            }


                                                        } //                                   
                                                    } else { // no RSC found so crossed pattern
                                                        problem_identified = true; // missingSequence between last_rightSC and last_leftSC    
                                                        sc_patterns.push_back(last_rightSC);
                                                        sc_patterns.push_back(last_leftSC);
                                                        //*********
                                                        // find overlap between LSC and RSC and make proper sequence to put here
                                                        // RSC + Overlap(if any) + LSC
                                                        //*********

                                                        if (!left_cons_seq.empty() && !right_cons_seq.empty()) {

                                                            /**16Dec21
                                                             filter missing sequence position using direction of mapped reads at that position*/

                                                            bool isMissingSequence = filter_missingSeq(contig_name, last_rightSC, bam_file);

                                                            if (isMissingSequence) {

                                                                string RSC_LSC = overlapSC(right_cons_seq, left_cons_seq);

                                                                if (!RSC_LSC.empty()) {
                                                                    // insert RSC_LSC at previous_last_rightSC or previous_last_leftSC
                                                                    out_file << contig_name << "\t" << "missingSequence (w overlap)" << "\t" << last_rightSC << "\t" << "forward" << "\t" << RSC_LSC << endl;
                                                                } else {
                                                                    //insert RSC + LSC at previous_last_rightSC or previous_last_leftSC
                                                                    out_file << contig_name << "\t" << "missingSequence (w/o overlap)" << "\t" << last_rightSC << "\t" << "forward" << "\t" << right_cons_seq << left_cons_seq << endl;
                                                                    }
                                                                crossed_found = true;
                                                                goto out_left_next;
                                                            } // missing sequence position filter
                                                        } // consensus sequence check
                                                    }
                                                    //  }
                                                }
                                            }
                                        }
                                        problem_identified = false;
                                        right_seq.clear();
                                        last_rightSC = 0;
                                        consective_rightSC = 0;
                                    }

                                }
                                problem_identified = false;
                                right_seq.clear();
                                last_rightSC = 0;
                                consective_rightSC = 0;
                            }

                            if (!problem_identified) { // else LSC seq not found it can be self/multigeneChimera or unsupported_insertion
                                int start, end;
                                first_frag = contig_seq.substr(0, last_leftSC);
                                
                                if (last_leftSC + 1 < contig_seq.size())
                                    sec_frag = contig_seq.substr(last_leftSC + 1, contig_seq.size() - 1); // fragment after
                                else
                                    sec_frag = "";
                                if (first_frag.size() < sec_frag.size()) { // take smallest fragment of contig to check for chimera or local mis-assembly
                                    query_frag = first_frag;
                                    first_frag = "";
                                    start = 0;
                                    end = last_leftSC;
                                } else { // sec fragment is smaller
                                    query_frag = sec_frag;
                                    sec_frag = "";
                                    start = last_leftSC;
                                    end = contig_seq.length() - 1;
                                }
                            
                                string target = first_frag + string(query_frag.size(), 'N') + sec_frag;
                                string blast_hit = perform_blast(query_frag, target, "");
                                vector<string> temp2;// = utils.split(blast_hit, "\t");
                                utils.str_split(blast_hit, temp2, delimiter);
                                blast_found = temp2[0];
                                int target_start = atoi(temp2[1].c_str());
                                int target_end = atoi(temp2[2].c_str());
                                bool repeat = false;
                              
                                if (!sc_patterns.empty()) { // check miss-assembled regions
                                    
                                    for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                       
                                        if (abs(target_start - *sc) < repeat_TH) {
                                           
                                            for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                
                                                if ((abs(target_end - *sc) < repeat_TH)) {
                                                    repeat = true;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }
                                if (repeat) {
                                    //                                                sc_patterns.push_back(target_start);
                                    //                                                sc_patterns.push_back(target_end);
                                    goto out_left_next;
                                } else if (blast_found != "*") {
                                    
                                    sc_patterns.push_back(target_start);
                                    sc_patterns.push_back(target_end);
                                    problem_identified = true; // selfChimera  start and end pos
                                    sc_patterns.push_back(start);
                                    sc_patterns.push_back(end);
                                    
                                    if (blast_found == "+") {
                                        out_file << contig_name << "\t" << "selfChimera" << "\t" << start << "-" << end << "\t" << "Forward" << "\t" << "Discard" << endl;
                                    } else {
                                        query_frag = utils.Rcomplement(query_frag);
                                        out_file << contig_name << "\t" << "selfChimera" << "\t" << start << "-" << end << "\t" << "\t" << "Reverse complement" << "Discard" << endl;
                                    }
                                    goto out_left_next;
                                } else {
                                   
                                    if (abs(end - start) >= ignore_short_seq) { // no cov check because its either one side of contig, not in the middle 
                                        
                                        if (!sc_patterns.empty()) { // check miss-assembled regions
                                           
                                            for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                              
                                                if (abs(start - *sc) < repeat_TH) {
                                                  
                                                    for (vector<int>::iterator sc = sc_patterns.begin(); sc != sc_patterns.end(); sc++) {
                                                        
                                                        if ((abs(end - *sc) < repeat_TH)) {
                                                            repeat = true;
                                                            break;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        if (repeat) {
                                            goto out_left_next;
                                        }
                                        //28oct21 avoid chimeras at the start or end positions of overlapped merged positions found in RI and merged by SCs
                                        float SC_TH_mc;

                                        if (find(right_SC_base_pos.begin(), right_SC_base_pos.end(), last_leftSC + 1) != right_SC_base_pos.end()) { // if crisscross pattern
                                            SC_TH_mc = SC_support_TH2;
                                        } else // single sided pattern
                                            SC_TH_mc = SC_support_TH;
                                        //cout <<float(consective_leftSC / cov.at(pos_check - 1)) << endl;
                                        if (float(consective_leftSC / cov.at(pos_check - 1)) * 100 >= SC_TH_mc) { // 10feb22
                                            if (merged_contigs_itr == merged_contigs.end()) {

                                                problem_identified = true;
                                                sc_patterns.push_back(start);
                                                sc_patterns.push_back(end);
                                                out_file << contig_name << "\t" << "multigeneChimera" << "\t" << start << "-" << end << "\t" << "Forward" << "\t" << "New contig" << endl;

                                            } else { //identified chimera positions must be at one side of merged points in the merged contigs

                                                if ((start < start_overlap_merged && end < start_overlap_merged) || (start > end_overlap_merged && end > end_overlap_merged)) {
                                                    problem_identified = true;
                                                    sc_patterns.push_back(start);
                                                    sc_patterns.push_back(end);
                                                    out_file << contig_name << "\t" << "multigeneChimera" << "\t" << start << "-" << end << "\t" << "Forward" << "\t" << "New contig" << endl;
                                                } // else go out
                                            }
                                        }
                                        goto out_left_next;

                                    }
                                }
                            }
                            // }
                        } //  consensus check
                    }// cov check
                }

out_left_next:
                ;
                consective_leftSC = 0;
                problem_identified = false;
                left_seq.clear();
                last_leftSC = 0;
            }
        }
    }

    //  }
    // cout<<"leaving SC"<<endl;
    return sc_patterns;
}

string mis_assembly_chimera::perform_blast(string query_seq, string subject1_seq, string subject2_seq) { // take two subject sequences
    utils utils;
    string entry;
    vector<string> temp;
    bool blast_found = false;
    int target_start, target_end;
    string strand = "*"; //  means no hit found
    string subject_seq = path_inter + "/subject_contig.fasta";
    string query_sequence = path_inter + "/query_seq.fasta";
    string blast_output = path_inter + "/blast_output";
    ofstream query, subject;
    query.open(query_sequence.c_str());
    subject.open(subject_seq.c_str());

    query << ">query" << endl << query_seq << endl;
    subject << ">subject1" << endl << subject1_seq << endl << ">subject2" << endl << subject2_seq << endl;

    query.close();
    subject.close();
    
    //string blast_command_left = exe_path + "external_tools/blast -stepSize=5 -repMatch=2253 -minScore=10 -minIdentity=50 " + subject_seq + " " + query_sequence + " " + blast_output + "  > /dev/null 2>&1";
    //system(blast_command_left.c_str());
    
    
    string overlap_byBLAST = "blastn -evalue 5 -word_size 4 -outfmt '7 qseqid qlen sseqid slen qcovhsp nident length mismatch qstart qend sstart send sstrand' -subject " + subject_seq  + " -query " + query_sequence + " -out " + blast_output + "  > /dev/null 2>&1";
    std::system(overlap_byBLAST.c_str());
   
    
    //cout << "blast done" << endl;
    ifstream output(blast_output.c_str());
    while (!(getline(output, entry).eof())) {

        if (entry[0] == '#') {
            continue;
        } else {
            utils.str_split(entry, temp, delimiter);
            string id1 = temp[0];
            string id2 = temp[2];
            int target_size = atoi(temp[3].c_str());

            target_start = atoi(temp[10].c_str());
            target_end = atoi(temp[11].c_str());

            float score1 = atoi(temp[5].c_str());
            float final_score = float(score1 / query_seq.size() ) * 100; // atoi(temp[4].c_str()); // score1 is number of matches

            if (final_score >= blast_score_TH) { // extracted seq found in contig hence selfChimera
                blast_found = true;

                if (temp[12] == "plus")
                    strand = "+";
                else // "minus"
                     strand = "-"; 

                goto out;
            }
        }
    }
out:
;
stringstream entry_new;
entry_new << strand << "\t" << target_start << "\t" << target_end;
string blast_hits = entry_new.str();
entry_new.str(string());
return blast_hits;
}

bool mis_assembly_chimera::filter_missingSeq(string contig_name, int last_rightSC, string bam_file) {
    utils utils;

    bool is_missingSeq = false, pad;
    
    std::size_t found = bam_file.find_last_of("/\\");
    string PATH = bam_file.substr(0, found);
    string bam_contig = PATH + "/" + contig_name + ".bam";
    
    stringstream start_region, end_region;

    start_region << last_rightSC - 100;
    end_region << last_rightSC + 100;
    // extract region of that contig from +-100 of identified position
    string extract_readData = "samtools view -b  " + bam_file + " ''" + contig_name + ":" + start_region.str() + "-" + end_region.str() + "'' " + " -o " + bam_contig + "  > /dev/null 2>&1";
  //  cout << extract_readData << endl;
    system(extract_readData.c_str()); //read s which have unmapped mates

    BamReader reader;
    BamTools::BamAlignment al;
    vector<CigarOp> cigar;
    BamTools::RefData ref_contig;
    vector<int> clipsize, read_pos, gen_pos;
    map<int, string> umappedMate_reads, dist_mapped_reads, same_dir_readpair, oppo_dir_readpair;
    
    int right_mapping_count = 0, left_mapping_count = 0, opposite_mapped_count = 0;
    // Opens a bam file
    if (!reader.Open(bam_contig)) {
        cerr << "Could not open BAM file." << endl;
        exit(0);
    }

    vector<RefData> ref;
    ref = reader.GetReferenceData();
    
    while (reader.GetNextAlignment(al)) {

        if (al.RefID >= 0 && al.IsMapped() && al.IsMateMapped()) { // to avoid unaligned reads which gives -1 RefID
            //ref_contig = ref.at(al.RefID);
            //cur_ref_id = al.RefID;
            // cigar = al.CigarData;

            if (al.GetSoftClips(clipsize, read_pos, gen_pos, pad) == true) {

                if (clipsize.size() == 1) { // if one sided clip , go for further checking

                    if (clipsize[0] == read_pos[0]) { // LSCs

                        // read and mate mapped on the same contig, SC position is same (or close as identified last_rightSC), read is in reverse direction and insert size is positive
                        if (al.RefID == al.MateRefID && abs(gen_pos[0] - last_rightSC) < 5 && (al.IsReverseStrand()) && (al.MatePosition + al.Length ) < last_rightSC ) { //al.InsertSize < 0) { //al.MatePosition > al.Position) {

                            left_mapping_count++;
                            opposite_mapped_count++;
                        }
                    } else if (clipsize[0] != read_pos[0]) { // RSCs

                        // read and mate mapped on the same contig, SC position is same (or close as identified last_rightSC), read is in forward direction and insert size is negative
                        if (al.RefID == al.MateRefID && abs(gen_pos[0] - last_rightSC) < 5 && !(al.IsReverseStrand()) && al.MatePosition > last_rightSC) {  //al.InsertSize > 0) {      

                            right_mapping_count++;
                            opposite_mapped_count++;
                        }
                    }
                } else {// if (clipsize.size() > 1) 

                    if (clipsize[0] == read_pos[0] && clipsize[1] <= max_SCs_oppositeSide_MCprocess) {

                        // read and mate mapped on the same contig, SC position is same (or close as identified last_rightSC), read is in reverse direction and insert size is positive
                        if (al.RefID == al.MateRefID && abs(gen_pos[0] - last_rightSC) < 5 && (al.IsReverseStrand()) && (al.MatePosition + al.Length ) < last_rightSC ) { //al.InsertSize < 0) { //al.MatePosition > al.Position) {

                            left_mapping_count++;
                            opposite_mapped_count++;
                        }
                    }
                    if (clipsize[1] != read_pos[1] && clipsize[0] <= max_SCs_oppositeSide_MCprocess) { //right SC clipsize[1] >= min_sc_len &&  25dec2020

                        // read and mate mapped on the same contig, SC position is same (or close as identified last_rightSC), read is in forward direction and insert size is negative
                        if (al.RefID == al.MateRefID && abs(gen_pos[1] - last_rightSC) < 5 && !(al.IsReverseStrand()) && al.MatePosition > last_rightSC) { //al.MatePosition > al.Position) {       

                            right_mapping_count++;
                            opposite_mapped_count++;
                        }
                    }
                }
                clipsize.clear();
                read_pos.clear();
                gen_pos.clear();
            }

        }
    }
    utils.remove_file(bam_contig);
    
    if(opposite_mapped_count >= min_sc_reads){//if (right_mapping_count >= min_sc_reads || left_mapping_count >= min_sc_reads ) { // min_sc_reads = 3
        is_missingSeq = true;
    }

    return is_missingSeq;
}

string mis_assembly_chimera::perform_blast(string query_seq, string subject1_seq) { // take single subject sequence
    utils utils;
    string entry;
    vector<string> temp;
    bool blast_found = false;
    int match_pos = -1; // can't initialize 0 b/c 0 can be real match position
    int end_pos;
    string subject_seq = path_inter + "/subject_contig.fasta";
    string query_sequence = path_inter + "/query_seq.fasta";
    string blast_output = path_inter + "/blast_output";
    ofstream query, subject;
    string strand = "*";
    query.open(query_sequence.c_str());
    subject.open(subject_seq.c_str());

    query << ">query" << endl << query_seq << endl;
    subject << ">subject" << endl << subject1_seq << endl;

    query.close();
    subject.close();
    
    //string blast_command_left = exe_path + "external_tools/blast -stepSize=2 -repMatch=2253 -minScore=10 -minIdentity=50 " + subject_seq + " " + query_sequence + " " + blast_output + "  > /dev/null 2>&1";
    //system(blast_command_left.c_str());
    
    string overlap_byBLAST = "blastn -evalue 5 -word_size 4 -outfmt '7 qseqid qlen sseqid slen qcovhsp nident length mismatch qstart qend sstart send sstrand' -subject " + subject_seq  + " -query " + query_sequence + " -out " + blast_output + "  > /dev/null 2>&1";
    std::system(overlap_byBLAST.c_str());

    //cout << "blast done" << endl;
    ifstream output(blast_output.c_str());
    while (!(getline(output, entry).eof())) {
        if (entry[0] == '#') {
            continue;
        } else {
            utils.str_split(entry, temp, delimiter);
            string id1 = temp[0];
            string id2 = temp[2];
            int target_size = atoi(temp[3].c_str());

            int target_start = atoi(temp[10].c_str());
            int target_end = atoi(temp[11].c_str());

            float score1 = atoi(temp[5].c_str());
            int final_score = float(score1 / query_seq.size()) * 100;//atoi(temp[4].c_str()); // score1 is number of matches

            if (final_score >= blast_score_TH) { // extracted seq found in contig hence selfChimera
                match_pos = target_start;
                end_pos = target_end;

                if (temp[12] == "plus")
                    strand = "+";
                else // minus
                    strand = "-"; 
                
                blast_found = true;
                goto out;
            }
        }
    }
out:
    ;
    stringstream entry_new;
    entry_new << match_pos << "\t" << strand << "\t" << end_pos;
    string blast_hits = entry_new.str();
    entry_new.str(string());
    output.close();
    return blast_hits;
}

string mis_assembly_chimera::overlapSC(string RSC, string LSC) {
    utils utils;
    int min_overlap_TH = 3;
    int pos_MI, pos_RI;
    string overlap = "", final, MI, RI, truncated_RI, truncated_MI;
    int overlapSize = 0;
    int i = 0, j = 0;
    char base_str1, base_str2;
    while (i < RSC.size()) {
        j = 0;
        while (j < LSC.size()) {
            if (RSC[i] == LSC[j]) {
                overlap = overlap + RSC[i];
                pos_MI = j; // size of overlap
                pos_RI = i;
                i++;
                j++;
            } else if (i == RSC.size()) // RSC processed 
                break;
            else {
                overlap = "";
                break;
            }

        }
        i++;
    }
    if (!(overlap.empty())) { //when overlap found till end
        RI = RSC.substr(0, pos_RI - overlap.size() + 1);
        MI = LSC.substr(pos_MI, LSC.size());
        truncated_RI = RSC.substr(pos_RI + 1, RSC.size() - 1); //end of contig1 left from overlap
        //std::transform(overlap.begin(), overlap.end(), overlap.begin(), ::tolower);
        if ((overlap.size() < min_overlap_TH) || (truncated_RI.size() > 1))
            final = "";
        else {
            final = RI + overlap + MI;
        }
    } else {
        //cout << "no overlap found" << endl;
    }
    return final;
}

int mis_assembly_chimera::filterPos_processFasta(string fastafile, string processed_fastafile, string identified_posFile, string mis_assembly_chimera_log, map<string, list<int> > &exon_boundary) {

    // cout << "in filter pos process Fasta" << endl;
    utils utils;
    std::size_t found; // = identified_posFile.find_last_of("/\\");
    //    string processed_fastafile = fastafile;
    // string final_pos = identified_posFile.substr(0, found) + "/final_considered_miassemblies4_1.fasta";

    ofstream processed_file;
    processed_file.open(processed_fastafile.c_str());

    ofstream final_pos_file;
    final_pos_file.open(mis_assembly_chimera_log.c_str());

    vector<int>::iterator st, en;
    vector<string> missing_seq, temp;
    missing_seq.reserve(200);
    
    vector<string> same_contig, processed_contigs;
    same_contig.reserve(10); processed_contigs.reserve(contig_count);
    
    list<string> ignore_cap3;

    string new_fasta_seq, fasta_seq, ID1 = "", ID2, pos_temp, ID_main;
    int start_pos, start_pos2, start_pos3, end_pos, end_pos2, end_pos3, insert_at, insert_at2, insert_at3, pos, count = 1, i = 0, split_count = 0, update_count = 0;
    
    int push = 0, pull = 0, push_pos = 0, push_till = 0, pull_till = 0, t_c_start, t_c_end;
    bool transchimera = false, translocation = false, inv_tra = false, ignore = false, self_chimera = false, to_push = true, at_exon_boundary = false, update = false;
    
    int contig_edge = 10; 
    
    list<int> ex_bo;
    map<string, list<int> >::iterator ex_bo_it;
    
    string entry, entry1;
    ifstream pot_chimera_file(identified_posFile.c_str()); // ID    type    start-end   strand  insert_pos/seq  Note

    vector <int> all_pos;
    vector <int> problem_type;
    while (!(getline(pot_chimera_file, entry).eof())) {

        utils.str_split(entry, temp, delimiter);
        if (count == 1) {
            ID1 = temp[0];
            same_contig.push_back(entry);
            
            if (found = temp[1].find("missingSequence") != string::npos)
                    missing_seq.push_back(temp[4]);
            //entry1 = entry;        
            count += 1;
            // update = false;
        } else { //second entry
            ID2 = temp[0];
            if (ID2 == ID1) {
                
                if (found = temp[1].find("missingSequence") != string::npos)
                    
                    missing_seq.push_back(temp[4]);

                same_contig.push_back(entry);
                entry1 = entry;
                ID1 = ID2;
                //update = false;

            } else { //entered in different contig, process previous
              //  cout << ID1 << endl;

                    
                
                map<string, string>::iterator it_fasta;
                it_fasta = AllFasta_data.find(ID1);
                
                 ex_bo_it = exon_boundary.find(ID1);
                    if (ex_bo_it != exon_boundary.end()) {
                        ex_bo = ex_bo_it->second;
                    }
                    
                fasta_seq = it_fasta->second; //utils.extract_fasta(ID1, fastafile);
                new_fasta_seq = fasta_seq;
                //    final_pos_file << ID1 << endl << new_fasta_seq << endl;

                get_problemData(same_contig, all_pos, problem_type);
                vector<int>::iterator it_pos = all_pos.begin();

                for (vector<int>::iterator it = problem_type.begin(); it != problem_type.end(); it++) {
                    if(!update){
                    if (*it == PROBLEM_MISSING_SEQ) {
                       // cout << " PROBLEM_MISSING_SEQ" << endl;
                        for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 

                            if (*it_pos != 0 && *it_pos >= *boundary - 5 && (*it_pos <= *boundary + 5)) {
                                at_exon_boundary = true;
                                it_pos++;
                                break;
                            }
                        }

                        if (!at_exon_boundary) {
                            string seq_to_insert = missing_seq.at(i);
                            i++;

                            // missing_seq.erase(missing_seq.begin() + 0);
                            push = seq_to_insert.size() - 1;
                            push_pos = *it_pos;
                            push_till = *max_element(all_pos.begin(), all_pos.end());
                            new_fasta_seq = new_fasta_seq.substr(0, push_pos - 1) + seq_to_insert + new_fasta_seq.substr(push_pos, new_fasta_seq.size() - 1);
                            final_pos_file << ID1 << " " << "m_s" << " " << *it_pos << endl;
                            // final_pos_file << new_fasta_seq << endl;
                            update_positions(all_pos, push_pos, push_till, push, to_push); // vector start_pos end_pos change_size push/pull
                            it_pos++;

                            update = true;
                            // cout << " PROBLEM_MISSING_SEQ.....end" << endl;
                        }
                    } else if (*it == PROBLEM_MIS_ASSEMBLY) {
                        //  cout << " PROBLEM_MIS_ASSEMBLY" << endl;
                        //usually at start or end of contig
                        // if start -> push by size of it
                        start_pos = *it_pos; // because reverse iterator
                        it_pos++;
                        end_pos = *it_pos;
                        it_pos++;
                        
                        for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 

                                if ((start_pos >= contig_edge && start_pos >= *boundary - 5 && (start_pos <= *boundary + 5)) || ( end_pos  <= fasta_seq.size() - contig_edge && end_pos >= *boundary - 5 && (end_pos <= *boundary + 5))) {
                                    goto next_problem;
                                }
                            }
                        
                        if(abs(start_pos - end_pos) < 10 ){ // 16 june 2021 to avaoid sc mc i and t of less than ten bases
                            continue;
                        }
                        pull = end_pos - start_pos;
                        //bool compare_positions(vector<int> all_positions, vector<string> problems, int p_start, int p_end, string problem);
                        ignore = compare_positions(all_pos, problem_type, start_pos, end_pos, *it);
                        if (!ignore) {
                            pull_till = *max_element(all_pos.begin(), all_pos.end());
                            if (start_pos == 0) { // miss assembly in the start of contig
                                new_fasta_seq = new_fasta_seq.substr(end_pos, new_fasta_seq.size() - 1);
                                update_positions(all_pos, end_pos, pull_till, pull, !to_push); // pull, vector start_pos end_pos change_size push/pull
                            } else if (end_pos == fasta_seq.size()) { // miss-assembly in the end of contig with push size
                                new_fasta_seq = new_fasta_seq.substr(0, start_pos  - 1);
                            } else {
                                new_fasta_seq = new_fasta_seq.substr(0, start_pos - 1) + new_fasta_seq.substr(end_pos, new_fasta_seq.size() - 1);
                                update_positions(all_pos, end_pos, pull_till, pull, !to_push); // pull, vector start_pos end_pos change_size push/pull
                            }

                            final_pos_file << ID1 << " " << "m_a" << " " << start_pos << " " << end_pos << endl;
                            update = true;
                         //   final_pos_file << new_fasta_seq << endl;
                        }
                        ignore = false;
                     //   cout << " PROBLEM_MIS_ASSEMBLY ...... end" <<endl;
                    }
                    else if (*it == PROBLEM_UNSUPPORTED_INSERTION) {
                      //  cout << " PROBLEM_UNSUPPORTED_INSERTION" <<endl;
                        start_pos = *it_pos;
                        it_pos++;
                        end_pos = *it_pos;
                        it_pos++;
                         if(abs(start_pos - end_pos) < 10 ){ // 16 june 2021
                            continue;
                        }
                              for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 

                                if ((start_pos >= contig_edge && start_pos >= *boundary - 5 && (start_pos <= *boundary + 5)) ||  ( end_pos  <= fasta_seq.size() - contig_edge && end_pos >= *boundary - 5 && (end_pos <= *boundary + 5))) {
                                    goto next_problem;
                                }
                            }
                        
                        
                        ignore = compare_positions(all_pos, problem_type, start_pos, end_pos, *it);
                        
                        if (!ignore) {
                            pull = end_pos - start_pos;
                            pull_till = *max_element(all_pos.begin(), all_pos.end());
                            new_fasta_seq = new_fasta_seq.substr(0, start_pos - 1) + new_fasta_seq.substr(end_pos - 1, new_fasta_seq.size() - 1);
                            update_positions(all_pos, end_pos, pull_till, pull, !to_push); // pull, vector start_pos end_pos change_size push/pull
                         //   final_pos_file << new_fasta_seq << endl;
                          final_pos_file << ID1 << " " << "u_i" << " " << start_pos << " " << end_pos << endl;
                          update = true;
                        }
                        ignore = false;
                     //   cout << " PROBLEM_UNSUPPORTED_INSERTION ...... end" <<endl;
                    }
                    else if (*it == PROBLEM_SELF_CHIMERA) {
                       // cout << " PROBLEM_SELF_CHIMERA" <<endl;
                        //1 even
                        //pull after end
                        start_pos = *it_pos;
                        it_pos++;
                        end_pos = *it_pos;
                        it_pos++;
                        
                         if(abs(start_pos - end_pos) < 10){
                            continue;
                        }
                              for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 

                                if ((start_pos >= contig_edge && start_pos >= *boundary - 5 && (start_pos <= *boundary + 5)) || ( end_pos  <= fasta_seq.size() - contig_edge && end_pos >= *boundary - 5 && (end_pos <= *boundary + 5))) { // && 12 dec21
                                    goto next_problem;
                                }
                            }
                        
                        
                        pull = end_pos - start_pos;
                        if (!self_chimera ) {
                            ignore = compare_positions(all_pos, problem_type, start_pos, end_pos, *it);
                            
                            if (!ignore) {
                                pull_till = *max_element(all_pos.begin(), all_pos.end());

                                if (start_pos == 0) { // self chimera in the start of contig
                                    new_fasta_seq = new_fasta_seq.substr(end_pos, new_fasta_seq.size() - 1);
                                    update_positions(all_pos, end_pos, pull_till, pull, !to_push); // pull,  vector start_pos end_pos change_size push/pull

                                    list<string>::iterator it_cleft;
                                    it_cleft = find(cap3_left_list.begin(), cap3_left_list.end(), ID1);

                                    /*Remove cap3 derivatives too if respective section get removed in self chimera*/
                                    if (it_cleft != cap3_left_list.end()) {
                                        string ID_cap3_left = ID1 + "_" + ROASTcap3_left_tag;
                                        ignore_cap3.push_back(ID_cap3_left);
                                    }

                                } else if (end_pos == fasta_seq.size()) { // self chimera in the end of contig with push size
                                    new_fasta_seq = new_fasta_seq.substr(0, start_pos - 1);

                                    list<string>::iterator it_cright;
                                    it_cright = find(cap3_right_list.begin(), cap3_right_list.end(), ID1);

                                    /*Remove cap3 derivatives too if respective section get removed in self chimera*/
                                    if (it_cright != cap3_right_list.end()) {
                                        string ID_cap3_right = ID1 + "_" + ROASTcap3_right_tag;
                                        ignore_cap3.push_back(ID_cap3_right);
                                    }

                                } else {
                                    new_fasta_seq = new_fasta_seq.substr(0, start_pos - 1) + new_fasta_seq.substr(end_pos - 1, new_fasta_seq.size() - 1);
                                    update_positions(all_pos, end_pos + 1, pull_till, pull, !to_push); // pull, vector start_pos end_pos change_size push/pull
                                }

                                // final_pos_file << new_fasta_seq << endl;
                                final_pos_file << ID1 << " " << "s_c" << " " << start_pos << " " << end_pos << endl;
                                update = true;
                                self_chimera = true;
                            }
                        }
                     //   cout << " PROBLEM_SELF_CHIMERA......end" <<endl;
                    } else if (*it == PROBLEM_MULTIGENE_CHIMERA) {
                      //  cout << " PROBLEM_MULTIGENE_CHIMERA " <<endl;
                        //1 even
                        //pull after end

                        if (!transchimera && (float(pull) / float(fasta_seq.size()) * 100 < 75)) {
                            // if (all_pos.size() > 2) {
                            start_pos = *it_pos;
                            t_c_start = std::distance(all_pos.begin(), it_pos);
                            it_pos++;
                            end_pos = *it_pos;
                            t_c_end = std::distance(all_pos.begin(), it_pos);
                            it_pos++;

                            if (abs(start_pos - end_pos) < 10) {
                                continue;
                            }                           
                            
                            for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 

                                if ((start_pos >= contig_edge && start_pos >= *boundary - 5 && (start_pos <= *boundary + 5)) || ( end_pos  <= fasta_seq.size() - contig_edge && end_pos >= *boundary - 5 && (end_pos <= *boundary + 5))) {
                                    goto next_problem;
                                }
                            }
                            
                            //ignore = compare_positions(all_pos, problem_type, start_pos, end_pos, *it);
                            
                            if (!ignore ) { //deal at the very end so new need to update other positions
                                pull = end_pos - start_pos;
                               // final_pos_file << ID1 << " " << "m_c" << " " << start_pos << " " << end_pos << endl;
                                transchimera = true;
                                update = true;
                            }
                            ignore = false;
                        } else {
                            start_pos = *it_pos; // because reverse iterator
                            it_pos++;
                            end_pos = *it_pos;
                            it_pos++;
                        }
                     //   cout << " PROBLEM_MULTIGENE_CHIMERA......end " <<endl;
                    } else if (*it == PROBLEM_TRANSLOCATION) {
                       // cout << " PPROBLEM_TRANSLOCATION " <<endl;
                        //1 odd
                        //pull values between end to insert

                        start_pos = *it_pos; // because reverse iterator
                        it_pos++;
                        end_pos = *it_pos;
                        it_pos++;
                        insert_at = *it_pos;
                        it_pos++;
                        
                         if(abs(start_pos - end_pos) < 10){
                            continue;
                        }
                        
                        for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 

                                if ((start_pos >= contig_edge && start_pos >= *boundary - 5 && (start_pos <= *boundary + 5)) || ( end_pos  <= fasta_seq.size() - contig_edge && end_pos >= *boundary - 5 && (end_pos <= *boundary + 5))) {
                                    goto next_problem;
                                }
                        }

                        pull = end_pos - start_pos;
                        int ind1 = 0;
                        if ((!translocation && !transchimera) && (float(pull) / float(fasta_seq.size()) * 100 < mis_assembly_chimera_par_len)) {
                           // ignore = compare_positions(all_pos, problem_type, start_pos, end_pos, insert_at, *it);

                            if (!ignore) {

                                if (insert_at < start_pos) {

                                    if (start_pos - insert_at <= ignore_short_seq) { //region between start position of translocation and insert size is less than 200 consider it as UI // 21Dec21
                                        
                                        end_pos = start_pos; // remove region between insert at and start position because of size
                                        start_pos = insert_at;
                                       
                                        pull = end_pos - start_pos;
                                        pull_till = *max_element(all_pos.begin(), all_pos.end());
                                        new_fasta_seq = new_fasta_seq.substr(0, insert_at) + new_fasta_seq.substr(start_pos, end_pos - start_pos) + new_fasta_seq.substr(end_pos, new_fasta_seq.size() - 1); // new_fasta_seq.substr(insert_at, start_pos - insert_at) + 
                                        update_positions(all_pos, end_pos, pull_till, pull, !to_push); // pull, vector start_pos end_pos change_size push/pull
                                        //update_positions(all_pos, insert_at, start_pos, pull, to_push);
                                        final_pos_file << ID1 << " " << "u_i" << " " << start_pos << " " << end_pos << endl;
                                        update = true;

                                    } else { // translocate normally

                                        new_fasta_seq = new_fasta_seq.substr(0, insert_at) + new_fasta_seq.substr(start_pos, end_pos - start_pos) + new_fasta_seq.substr(insert_at, start_pos - insert_at) + new_fasta_seq.substr(end_pos, new_fasta_seq.size() - 1);
                                        update_positions(all_pos, insert_at, start_pos, pull, to_push);
                                        final_pos_file << ID1 << " " << "t" << " " << start_pos << " " << end_pos << " " << insert_at << endl;
                                        translocation = true;
                                        update = true;
                                    }

                                } else if (insert_at > end_pos) {

                                    if (insert_at - end_pos <= ignore_short_seq) {
                                       
                                        start_pos = end_pos; // remove region between insert at and start position because of size
                                        end_pos = insert_at;
                                        
                                        pull = end_pos - start_pos;
                                        pull_till = *max_element(all_pos.begin(), all_pos.end());
                                        new_fasta_seq = new_fasta_seq.substr(0, start_pos) + new_fasta_seq.substr(start_pos, end_pos - start_pos) + new_fasta_seq.substr(insert_at, new_fasta_seq.size() - 1); //+ new_fasta_seq.substr(end_pos, insert_at - end_pos)
                                        update_positions(all_pos, end_pos, pull_till, pull, !to_push); // pull, vector start_pos end_pos change_size push/pull
                                        //update_positions(all_pos, end_pos, insert_at, pull, !to_push); //"pull");
                                        final_pos_file << ID1 << " " << "u_i" << " " << start_pos << " " << end_pos << endl;
                                        update = true;
                                        
                                    } else {
                                        
                                        new_fasta_seq = new_fasta_seq.substr(0, start_pos) + new_fasta_seq.substr(end_pos, insert_at - end_pos) + new_fasta_seq.substr(start_pos, end_pos - start_pos) + new_fasta_seq.substr(insert_at, new_fasta_seq.size() - 1);
                                        update_positions(all_pos, end_pos, insert_at, pull, !to_push); //"pull");
                                        final_pos_file << ID1 << " " << "t" << " " << start_pos << " " << end_pos << " " << insert_at << endl;
                                        translocation = true;
                                        update = true;
                                    }
                                }

                                
                                //   final_pos_file << new_fasta_seq << endl;
                            }
                        }
                        ignore = false;
                        //   cout << " PPROBLEM_TRANSLOCATION.....end " <<endl;
                    } else if (*it == PROBLEM_INVERSION) {
                        //  cout << " PROBLEM_INVERSION" <<endl;
                        //1 even
                        start_pos = *it_pos; // because reverse iterator
                        it_pos++;
                        end_pos = *it_pos;
                        it_pos++;
                        
                         if(abs(start_pos - end_pos) < 10){
                            continue;
                        }
                        
                        for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 

                                if ((start_pos >= contig_edge && start_pos >= *boundary - 5 && (start_pos <= *boundary + 5)) || ( end_pos  <= fasta_seq.size() - contig_edge && end_pos >= *boundary - 5 && (end_pos <= *boundary + 5))) {
                                    goto next_problem;
                                }
                            }
                        
                        pull = end_pos - start_pos;
                       // ignore = compare_positions(all_pos, problem_type, start_pos, end_pos, *it);
                        if (!ignore) {
                            string fasta_rc = new_fasta_seq.substr(start_pos, end_pos - start_pos);
                            fasta_rc = utils.Rcomplement(fasta_rc);
                            new_fasta_seq = new_fasta_seq.substr(0, start_pos) + fasta_rc + new_fasta_seq.substr(end_pos, new_fasta_seq.size() - 1);
                           // final_pos_file << new_fasta_seq << endl;
                            final_pos_file << ID1 << " " << "i" << " " << start_pos << " " << end_pos << endl;
                            update = true;
                        }
                        ignore = false;
                     //   cout << " PROBLEM_INVERSION.......end" <<endl;
                    } else if (*it == PROBLEM_INVERSION_TRANSLOCATION) {
                      //  cout << " PROBLEM_INVERSION_TRANSLOCATION" <<endl;
                        //1 odd
                        //pull values between end to insert
                        start_pos = *it_pos;
                        it_pos++;
                        end_pos = *it_pos;
                        it_pos++;
                        insert_at = *it_pos;
                        it_pos++;
                        
                        if(abs(start_pos - end_pos) < 10 ){
                            continue;
                            }

                            for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 

                                if ((start_pos >= contig_edge && start_pos >= *boundary - 5 && (start_pos <= *boundary + 5)) || (end_pos <= fasta_seq.size() - contig_edge && end_pos >= *boundary - 5 && (end_pos <= *boundary + 5)) && (insert_at >= contig_edge && insert_at <= fasta_seq.size() - contig_edge && insert_at >= *boundary - 5 && (insert_at <= *boundary + 5))) {
                                    goto next_problem;
                                }
                            }

                            pull = end_pos - start_pos;
                            int ind1 = 0;

                            if (!translocation && !transchimera && (float(pull) / float(fasta_seq.size()) * 100 < mis_assembly_chimera_par_len)) {

                                // ignore = compare_positions(all_pos, problem_type, start_pos, end_pos, insert_at, *it);

                                if (!ignore) {

                                    string fasta_rc = new_fasta_seq.substr(start_pos, end_pos - start_pos);
                                    fasta_rc = utils.Rcomplement(fasta_rc);

                                    if (insert_at < start_pos) {
                                        if (start_pos - insert_at <= ignore_short_seq) { //region between start position of translocation and insert size is less than 200 consider it as UI // 21Dec21

                                            end_pos = start_pos; // remove region between insert at and start position because of size
                                            start_pos = insert_at;

                                            pull = end_pos - start_pos;
                                            pull_till = *max_element(all_pos.begin(), all_pos.end());
                                            new_fasta_seq = new_fasta_seq.substr(0, insert_at) + new_fasta_seq.substr(start_pos, end_pos - start_pos) + new_fasta_seq.substr(end_pos, new_fasta_seq.size() - 1); // new_fasta_seq.substr(insert_at, start_pos - insert_at) + 
                                            update_positions(all_pos, end_pos, pull_till, pull, !to_push); // pull, vector start_pos end_pos change_size push/pull
                                            //update_positions(all_pos, insert_at, start_pos, pull, to_push);
                                            final_pos_file << ID1 << " " << "u_i" << " " << start_pos << " " << end_pos << endl;
                                            update = true;

                                        } else { // translocate normally

                                            new_fasta_seq = new_fasta_seq.substr(0, insert_at) + fasta_rc + new_fasta_seq.substr(insert_at, start_pos - insert_at) + new_fasta_seq.substr(end_pos, new_fasta_seq.size() - 1);
                                            update_positions(all_pos, insert_at, start_pos, pull, to_push);
                                            translocation = true;
                                            update = true;
                                            final_pos_file << ID1 << " " << "t_i" << " " << start_pos << " " << end_pos << " " << insert_at << endl;
                                        }
                                    } else if (insert_at > end_pos) {
                                        if (insert_at - end_pos <= ignore_short_seq) {

                                            start_pos = end_pos; // remove region between insert at and start position because of size
                                            end_pos = insert_at;

                                            pull = end_pos - start_pos;
                                            pull_till = *max_element(all_pos.begin(), all_pos.end());
                                            new_fasta_seq = new_fasta_seq.substr(0, start_pos) + new_fasta_seq.substr(start_pos, end_pos - start_pos) + new_fasta_seq.substr(insert_at, new_fasta_seq.size() - 1); //+ new_fasta_seq.substr(end_pos, insert_at - end_pos)
                                            update_positions(all_pos, end_pos, pull_till, pull, !to_push); // pull, vector start_pos end_pos change_size push/pull
                                            //update_positions(all_pos, end_pos, insert_at, pull, !to_push); //"pull");
                                            final_pos_file << ID1 << " " << "u_i" << " " << start_pos << " " << end_pos << endl;
                                            update = true;

                                        } else {
                                            new_fasta_seq = new_fasta_seq.substr(0, start_pos) + new_fasta_seq.substr(end_pos, insert_at - end_pos - 1) + fasta_rc + new_fasta_seq.substr(insert_at, new_fasta_seq.size() - 1);
                                            update_positions(all_pos, end_pos, insert_at, pull, !to_push); // "pull");
                                            translocation = true;
                                            update = true;
                                            final_pos_file << ID1 << " " << "t_i" << " " << start_pos << " " << end_pos << " " << insert_at << endl;
                                        }
                                    }
                                    //  final_pos_file << new_fasta_seq << endl;
                                }
                            }
                            ignore = false;
                            //  cout << " PROBLEM_INVERSION_TRANSLOCATION......end" <<endl;
                    } else if (*it == PROBLEM_INVERSION_TRANSLOCATION_TRANSLOCATION) {
                        //  cout << " PROBLEM_INVERSION_TRANSLOCATION_TRANSLOCATION" <<endl;
                        //2 even
                        //pull values between end1 to insert1 & end2 insert2
                        start_pos = *it_pos;
                        it_pos++;
                        end_pos = *it_pos;
                        it_pos++;
                        insert_at = *it_pos;
                        it_pos++;
                        start_pos2 = *it_pos;
                        it_pos++;
                        end_pos2 = *it_pos;
                        it_pos++;
                        insert_at2 = *it_pos;
                        it_pos++;
                        
                         if(abs(start_pos - end_pos) < 10 || abs(start_pos2 - end_pos2) < 10){
                            continue;
                        }
                        for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 

                                if ((start_pos >= contig_edge && start_pos >= *boundary - 5 && (start_pos <= *boundary + 5)) || ( end_pos  <= fasta_seq.size() - contig_edge && end_pos >= *boundary - 5 && (end_pos <= *boundary + 5)) && (insert_at >= contig_edge && insert_at <= fasta_seq.size() - contig_edge && insert_at >= *boundary - 5 && (insert_at <= *boundary + 5))) {
                                    goto next_problem;
                                }
                            }
                        for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 

                                if ((start_pos2 >= contig_edge && start_pos2 >= *boundary - 5 && (start_pos2 <= *boundary + 5)) || ( end_pos2  <= fasta_seq.size() - contig_edge && end_pos2 >= *boundary - 5 && (end_pos2 <= *boundary + 5)) && (insert_at2 >= contig_edge && insert_at2 <= fasta_seq.size() - contig_edge && insert_at2 >= *boundary - 5 && (insert_at2 <= *boundary + 5))) {
                                    goto next_problem;
                                }
                            }
                        
                        pull = end_pos - start_pos;
                        if (!translocation && !transchimera && (float(pull) / float(fasta_seq.size()) * 100 < mis_assembly_chimera_par_len)) {
                         //   ignore = compare_positions(all_pos, problem_type, start_pos, end_pos, insert_at, *it);
                               bool ignore2 = false; //compare_positions(all_pos, problem_type, start_pos2, end_pos2, insert_at2, *it);
                            
                            if (!ignore && !ignore2) {
                           
                                string fasta_rc = new_fasta_seq.substr(start_pos, end_pos - start_pos);
                                fasta_rc = utils.Rcomplement(fasta_rc);
                                
                                if (insert_at < start_pos) {
                                    new_fasta_seq = new_fasta_seq.substr(0, insert_at) + fasta_rc + new_fasta_seq.substr(insert_at, start_pos - insert_at) + new_fasta_seq.substr(end_pos, new_fasta_seq.size() - 1);
                                    update_positions(all_pos, insert_at, start_pos, pull, to_push);
                           
                                } else if (insert_at > end_pos) {
                                    new_fasta_seq = new_fasta_seq.substr(0, start_pos) + new_fasta_seq.substr(end_pos, insert_at - end_pos - 1) + fasta_rc + new_fasta_seq.substr(insert_at, new_fasta_seq.size() - 1);
                                    update_positions(all_pos, end_pos, insert_at, pull, !to_push); //pull
                                }
                                translocation = true;
                                update = true;
                              //  final_pos_file << new_fasta_seq << endl;
                                final_pos_file << ID1 << " " << "i_t_t" << " " << start_pos << " " << end_pos << " " << insert_at << endl;
                                //sec
                                int pull2 = end_pos2 - start_pos2;
                                int ind1 = 0;
                           
                                if (insert_at2 < start_pos2) {
                           
                                    new_fasta_seq = new_fasta_seq.substr(0, insert_at2) + new_fasta_seq.substr(start_pos2, end_pos2 - start_pos2) + new_fasta_seq.substr(insert_at2, start_pos2 - insert_at2) + new_fasta_seq.substr(end_pos2, new_fasta_seq.size() - 1);
                                    update_positions(all_pos, insert_at2, start_pos2, pull2, to_push);
                            
                                } else if (insert_at2 > end_pos2) {
                           
                                    new_fasta_seq = new_fasta_seq.substr(0, start_pos2) + new_fasta_seq.substr(end_pos2, insert_at2 - end_pos2 - 1) + new_fasta_seq.substr(start_pos2, end_pos2 - start_pos2) + new_fasta_seq.substr(insert_at2, new_fasta_seq.size() - 1);
                                    update_positions(all_pos, end_pos2, insert_at2, pull2, !to_push) ; //"pull");
                                }
                                translocation = true;

                              //  final_pos_file << new_fasta_seq << endl;
                                final_pos_file << ID1 << " " << "i_t_t"<< " " << start_pos2 << " " << end_pos2 << " " << insert_at2 << endl;
                            }
                        }
                      //  cout << " PROBLEM_INVERSION_TRANSLOCATION_TRANSLOCATION......end" <<endl;
                    } else if (*it == PROBLEM_INVERSION_TRANSLOCATION_SELFCHIMERA) {
                      //  cout << " PROBLEM_INVERSION_TRANSLOCATION_SELFCHIMERA" <<endl;
                        //2 odd
                        //pull values between end1 to insert1 
                        //pull after end2    by size of                  
                        start_pos = *it_pos;
                        it_pos++;
                        end_pos = *it_pos;
                        it_pos++;
                        insert_at = *it_pos;
                        it_pos++;
                        //second
                        start_pos2 = *it_pos;
                        it_pos++;
                        end_pos2 = *it_pos;
                        it_pos++;
                        pull = end_pos - start_pos;

                        if (abs(start_pos - end_pos) < 10 || abs(start_pos2 - end_pos2) < 10) {
                            continue;
                        }
                        for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 

                            if ((start_pos!= 0 && start_pos >= *boundary - 5 && (start_pos <= *boundary + 5)) || ( end_pos  <= fasta_seq.size() - contig_edge && end_pos >= *boundary - 5 && (end_pos <= *boundary + 5)) || (insert_at >= contig_edge && insert_at <=  fasta_seq.size() - contig_edge && insert_at >= *boundary - 5 && (insert_at <= *boundary + 5))) {
                                goto next_problem;
                            }
                        }
                        for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 

                            if ((start_pos2 >= contig_edge && start_pos2 >= *boundary - 5 && (start_pos2 <= *boundary + 5)) || ( end_pos2  <= fasta_seq.size() - contig_edge && end_pos2 >= *boundary - 5 && (end_pos2 <= *boundary + 5))) {
                                goto next_problem;
                            }
                        }

                        if (!translocation && !transchimera && (float(pull) / float(fasta_seq.size()) * 100 < mis_assembly_chimera_par_len)) {
                        //    ignore = compare_positions(all_pos, problem_type, start_pos, end_pos, insert_at, *it);
                           bool ignore2 = false; //compare_positions(all_pos, problem_type, start_pos2, end_pos2, *it);
                            
                            if (!ignore && !ignore2) {
                                string fasta_rc = new_fasta_seq.substr(start_pos, end_pos - start_pos);
                                fasta_rc = utils.Rcomplement(fasta_rc);

                                if (insert_at < start_pos) {
                                    new_fasta_seq = new_fasta_seq.substr(0, insert_at) + fasta_rc + new_fasta_seq.substr(insert_at, start_pos - insert_at) + new_fasta_seq.substr(end_pos, new_fasta_seq.size() - 1);
                                    update_positions(all_pos, insert_at, start_pos, pull, to_push);
                                } else if (insert_at > end_pos) {
                                    new_fasta_seq = new_fasta_seq.substr(0, start_pos) + new_fasta_seq.substr(end_pos, insert_at - end_pos - 1) + fasta_rc + new_fasta_seq.substr(insert_at, new_fasta_seq.size() - 1);
                                    update_positions(all_pos, end_pos, insert_at, pull, !to_push) ;//"pull");
                                }
                                translocation = true;
                                update = true;
                               // final_pos_file << new_fasta_seq << endl;
                                final_pos_file << ID1 << " " << "i_t_sc" << " " << start_pos << " " << end_pos << " " << insert_at <<  endl;
                                // new_fasta_seq = new_fasta_seq.substr(0, start_pos - 1) + new_fasta_seq.substr(end_pos + 1, insert_at - end_pos - 1) + fasta_rc + new_fasta_seq.substr(insert_at, fasta_seq.size());
                                int pull2 = end_pos2 - start_pos2;
                                pull_till = *max_element(all_pos.begin(), all_pos.end());
                                new_fasta_seq = new_fasta_seq.substr(0, start_pos2 - 1) + new_fasta_seq.substr(end_pos2 - 1, new_fasta_seq.size() - 1);
                                update_positions(all_pos, end_pos, pull_till, pull2, !to_push) ;//"pull"); // vector start_pos end_pos change_size push/pull
                                //final_pos_file << new_fasta_seq << endl;
                                final_pos_file << ID1 << " " << "i_t_sc" << " " << start_pos2 << " " << end_pos2 << " " << insert_at2 << endl;
                            }
                        }
                        ignore = false;
                        
                       // cout << " PROBLEM_INVERSION_TRANSLOCATION_SELFCHIMERA........end" <<endl;
                    } else if (*it == PROBLEM_INVERSION_TRANSLOCATION_INVERSION_TRANSLOCATION) {
                      //  cout << " PROBLEM_INVERSION_TRANSLOCATION_INVERSION_TRANSLOCATION" <<endl;
                        
                        //2 even
                        //pull values between end1 to insert1 & end2 insert2
                        start_pos = *it_pos; // because reverse iterator
                        it_pos++;
                        end_pos = *it_pos;
                        it_pos++;
                        insert_at = *it_pos;
                        it_pos++;
                        pull = end_pos - start_pos;
                        // second
                        start_pos2 = *it_pos; // because reverse iterator
                        it_pos++;
                        end_pos2 = *it_pos;
                        it_pos++;
                        insert_at2 = *it_pos;
                        it_pos++;
                        int pull2 = end_pos2 - start_pos2;
                        int ind1 = 0;

                        if (abs(start_pos - end_pos) < 10 || abs(start_pos2 - end_pos2) < 10) {
                            continue;
                        }
                        for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 

                            if ((start_pos >= contig_edge && start_pos >= *boundary - 5 && (start_pos <= *boundary + 5)) || (end_pos <= fasta_seq.size() - contig_edge && end_pos >= *boundary - 5 && (end_pos <= *boundary + 5)) || (insert_at >= contig_edge && insert_at <= fasta_seq.size() - contig_edge && insert_at >= *boundary - 5 && (insert_at <= *boundary + 5))) {
                                goto next_problem;
                            }
                        }
                        for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 

                            if ((start_pos2 >= contig_edge && start_pos2 >= *boundary - 5 && (start_pos2 <= *boundary + 5)) || (end_pos2 <= fasta_seq.size() - contig_edge && end_pos2 >= *boundary - 5 && (end_pos2 <= *boundary + 5)) || (insert_at2 >= contig_edge && insert_at2 <= fasta_seq.size() - contig_edge && insert_at2 >= *boundary - 5 && (insert_at2 <= *boundary + 5))) {
                                goto next_problem;
                            }
                        }

                        if (!translocation && !transchimera && (float(pull) / float(fasta_seq.size()) * 100 < mis_assembly_chimera_par_len)) {

                            //                            vector<int>::iterator test = all_pos.end();
                            //                                --test;
                           // ignore = compare_positions(all_pos, problem_type, start_pos, end_pos, insert_at, *it);
                            bool ignore2 = true;//compare_positions(all_pos, problem_type, start_pos2, end_pos2, insert_at2, *it);
                            
                            if (!ignore && !ignore2) {
                                string fasta_rc = new_fasta_seq.substr(start_pos, end_pos - start_pos);
                                fasta_rc = utils.Rcomplement(fasta_rc);
                                
                                if (insert_at < start_pos) {
                                    new_fasta_seq = new_fasta_seq.substr(0, insert_at) + fasta_rc + new_fasta_seq.substr(insert_at, start_pos - insert_at) + new_fasta_seq.substr(end_pos, new_fasta_seq.size() - 1);
                                    //    all_pos = update_positions(all_pos, insert_at, start_pos, pull, "push");
                                } else if (insert_at > end_pos) {
                                    new_fasta_seq = new_fasta_seq.substr(0, start_pos) + new_fasta_seq.substr(end_pos, insert_at - end_pos - 1) + fasta_rc + new_fasta_seq.substr(insert_at, new_fasta_seq.size() - 1);
                                    //   all_pos = update_positions(all_pos, end_pos, insert_at, pull, "pull");
                                }
                                translocation = true;
                                update = true;
                                //final_pos_file << new_fasta_seq << endl;
                                final_pos_file << ID1 << " " << "i_t_i_t" << " " << start_pos << " " << end_pos << " " << insert_at << endl;

                                fasta_rc = new_fasta_seq.substr(start_pos2, end_pos2 - start_pos2);
                                fasta_rc = utils.Rcomplement(fasta_rc);
                                ind1 = 0;
                                if (insert_at2 < start_pos2) {
                                    new_fasta_seq = new_fasta_seq.substr(0, insert_at) + fasta_rc + new_fasta_seq.substr(insert_at, start_pos - insert_at) + new_fasta_seq.substr(end_pos, new_fasta_seq.size() - 1);
                                    update_positions(all_pos, insert_at2, start_pos2, pull2, to_push);
                                } else if (insert_at2 > end_pos2) {
                                    new_fasta_seq = new_fasta_seq.substr(0, start_pos) + new_fasta_seq.substr(end_pos, insert_at - end_pos - 1) + fasta_rc + new_fasta_seq.substr(insert_at, new_fasta_seq.size() - 1);
                                    update_positions(all_pos, end_pos2, insert_at2, pull2, !to_push) ;//"pull");
                                }
                                translocation = true;
                                update_positions(all_pos, end_pos, insert_at, pull, !to_push) ;// "pull");
                               // final_pos_file << new_fasta_seq << endl;
                                final_pos_file << ID1 << " " << *it << "i_t_i_t" << " " << start_pos2 << " " << end_pos2 << " " << insert_at2 << endl;
                            }
                        }
                      //  cout << " PROBLEM_INVERSION_TRANSLOCATION_INVERSION_TRANSLOCATION...........end" <<endl;
                        
                    } else if (*it == PROBLEM_INVERSION_TRANSLOCATION_MULTIGENECHIMERA) {
                        
                      //  cout << " PROBLEM_INVERSION_TRANSLOCATION_MULTIGENECHIMERA" <<endl;
                        //2 odd
                        //pull values between end1 to insert1 
                        //pull after end2 
                        start_pos = *it_pos;
                        it_pos++;
                        end_pos = *it_pos;
                        it_pos++;
                        insert_at = *it_pos;
                        it_pos++;
                        pull = end_pos - start_pos;

                        if (abs(start_pos - end_pos) < 10) {
                            continue;
                        }
                        for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 

                            if ((start_pos >= contig_edge && start_pos >= *boundary - 5 && (start_pos <= *boundary + 5)) || ( end_pos  <= fasta_seq.size() - contig_edge && end_pos >= *boundary - 5 && (end_pos <= *boundary + 5)) || (insert_at >= contig_edge && insert_at <= fasta_seq.size() - contig_edge && insert_at >= *boundary - 5 && (insert_at <= *boundary + 5))) {
                                goto next_problem;
                            }
                        }

                        if (!translocation && !transchimera && (float(pull) / float(fasta_seq.size()) * 100 < mis_assembly_chimera_par_len)) {
                            string fasta_rc = new_fasta_seq.substr(start_pos, end_pos - start_pos);
                            fasta_rc = utils.Rcomplement(fasta_rc);
                            //ignore = compare_positions(all_pos, problem_type, start_pos, end_pos, insert_at, *it);
                            
                            if (!ignore) {
                                if (insert_at < start_pos) {
                                    new_fasta_seq = new_fasta_seq.substr(0, insert_at) + fasta_rc + new_fasta_seq.substr(insert_at, start_pos - insert_at) + new_fasta_seq.substr(end_pos, new_fasta_seq.size() - 1);
                                    update_positions(all_pos, insert_at, start_pos, pull, to_push);
                            
                                } else if (insert_at > end_pos) {
                                    new_fasta_seq = new_fasta_seq.substr(0, start_pos) + new_fasta_seq.substr(end_pos, insert_at - end_pos - 1) + fasta_rc + new_fasta_seq.substr(insert_at, new_fasta_seq.size() - 1);
                                    update_positions(all_pos, end_pos, insert_at, pull, !to_push) ;//"pull");
                                }
                                translocation = true;
                                update = true;
                              //  final_pos_file << new_fasta_seq << endl;
                                final_pos_file << ID1 << " " << "i_t_m" << " " << start_pos << " " << end_pos << " " << insert_at << endl;

                                //second transchimera
                                start_pos2 = *it_pos; // because reverse iterator
                                t_c_start = std::distance(all_pos.begin(), it_pos);
                                it_pos++;
                                end_pos2 = *it_pos;
                                t_c_end = std::distance(all_pos.begin(), it_pos);
                                it_pos++;
                                int pull2 = end_pos2 - start_pos2;
                            } else {
                                start_pos = *it_pos; // because reverse iterator
                                it_pos++;
                                end_pos = *it_pos;
                                it_pos++;
                            }
                        } else { // ignore and remove entry
                            start_pos = *it_pos; // because reverse iterator
                            it_pos++;
                            end_pos = *it_pos;
                            it_pos++;
                        }
                        
                     //   cout << " PROBLEM_INVERSION_TRANSLOCATION_MULTIGENECHIMERA..........end" <<endl;
                    } else if (*it == PROBLEM_INVERSION_TRANSLOCATIONandSELFCHIMERA) {
                        
                      //  cout << " PROBLEM_INVERSION_TRANSLOCATIONandSELFCHIMERA" <<endl;
                        //3 even
                        //leave first as it is as nothing will come between facing pattern which inverse and translocate
                        //pull values between end2 and insert2
                        //insert_at == start2
                        start_pos = *it_pos;
                        it_pos++;
                        end_pos = *it_pos;
                        it_pos++;
                        insert_at = *it_pos;
                        it_pos++;
                        pull = end_pos - start_pos;
                        start_pos2 = *it_pos;
                        it_pos++;
                        end_pos2 = *it_pos;
                        it_pos++;
                        insert_at2 = *it_pos;
                        it_pos++;
                        int pull2 = end_pos2 - start_pos2;

                        if (abs(start_pos - end_pos) < 10 || abs(start_pos2 - end_pos2) < 10) {
                            continue;
                        }
                        for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 

                            if ((start_pos >= contig_edge && start_pos >= *boundary - 5 && (start_pos <= *boundary + 5)) || ( end_pos  <= fasta_seq.size() - contig_edge && end_pos >= *boundary - 5 && (end_pos <= *boundary + 5)) || (insert_at >= contig_edge && insert_at <=  fasta_seq.size() - contig_edge && insert_at >= *boundary - 5 && (insert_at <= *boundary + 5))) {
                                goto next_problem;
                            }
                        }
                        for (list<int>::iterator boundary = ex_bo.begin(); boundary != ex_bo.end(); boundary++) { // to avoid selecting soft clips around exon-exon boundaries 

                            if ((start_pos2 >= contig_edge && start_pos2 >= *boundary - 5 && (start_pos2 <= *boundary + 5)) || ( end_pos2  <= fasta_seq.size() - contig_edge && end_pos2 >= *boundary - 5 && (end_pos2 <= *boundary + 5)) || (insert_at2 >= contig_edge && insert_at2 <=  fasta_seq.size() - contig_edge && insert_at2 >= *boundary - 5 && (insert_at2 <= *boundary + 5))) {
                                goto next_problem;
                            }
                        }
                        
                        if (!translocation && !transchimera && (float(pull) / float(fasta_seq.size()) * 100 < mis_assembly_chimera_par_len)) {
                            
                           // ignore = compare_positions(all_pos, problem_type, start_pos, end_pos, insert_at, *it);
                            bool ignore2 = true; // compare_positions(all_pos, problem_type, start_pos2, end_pos2, insert_at2, *it);
                            
                            if (!ignore && !ignore2) {
                                
                                string fasta_rc1 = new_fasta_seq.substr(start_pos, end_pos - start_pos);
                                fasta_rc1 = utils.Rcomplement(fasta_rc1);
                                string fasta_rc2 = new_fasta_seq.substr(start_pos2, end_pos2 - start_pos2);
                                fasta_rc2 = utils.Rcomplement(fasta_rc2);
                                
                                if (insert_at < start_pos) {
                                    update_positions(all_pos, insert_at, start_pos, pull, to_push);
                                } else if (insert_at > end_pos) {
                                    update_positions(all_pos, end_pos, insert_at, pull, !to_push) ;//"pull");
                                }
                                // for secnd inversion translocation
                                if (insert_at2 < start_pos2) {
                                    update_positions(all_pos, insert_at, start_pos, pull, to_push);
                                }
                                else if (insert_at2 > end_pos2) {
                                    update_positions(all_pos, end_pos, insert_at, pull, !to_push) ;//"pull");
                                }
                                new_fasta_seq = fasta_rc2 + fasta_rc1; //remove middle part simply

                                // for third selfchimera
                                start_pos3 = *it_pos; // because reverse iterator
                                it_pos++;
                                end_pos3 = *it_pos;
                                it_pos++;
                                int pull3 = end_pos3 - start_pos3;
                                pull_till = *max_element(all_pos.begin(), all_pos.end());
                                int ind = 0;
                                update_positions(all_pos, end_pos3, pull_till - 1, pull3, !to_push); //"pull");
                                translocation = true;
                                update = true;
                                it_pos++;
                                //    final_pos_file << new_fasta_seq << endl;
                                final_pos_file << ID1 << " " << "i_t&sc" << " " << start_pos << " " << end_pos << " " << insert_at << endl;

                                //         cout << " PROBLEM_INVERSION_TRANSLOCATIONandSELFCHIMERA...............end" <<endl;
                            }
                        }
                    }
next_problem:
                    ;
                    at_exon_boundary = false;
                    pull = 0;
                }
            }
                ex_bo.clear();
                if (update) { // not filtered by exon-exon boundary check
                    size_t update_ind = ID1.find(tool_name_flag);
                    size_t update_ind2 = ID1.find_last_of(tool_name_flag);
                    string header1, header2;

                    ID_main = ID1;
                    std::list<std::string>::iterator it;
                    it = std::find(split_contigs.begin(), split_contigs.end(), ID_main);

                    if (transchimera && it == split_contigs.end()) {//  &&  ID1.find("a00") == string::npos && ID1.find("b00") == string::npos ){
                        //   cout << "start splitting multigene chimera" << endl;
                        string new_contig, new_contog2;

                        if (all_pos[t_c_start] == 0) //transchimera in the start of contig
                        {
                            new_contig = new_fasta_seq.substr(all_pos[t_c_start], all_pos[t_c_end]);
                            new_contog2 = new_fasta_seq.substr(all_pos[t_c_end], new_fasta_seq.size());

                        } else if (all_pos[t_c_end] >= new_fasta_seq.size() - 5) //transchimera in the end of contig // -5 added on 11 May due to bug
                        {
                            new_contig = new_fasta_seq.substr(all_pos[t_c_start], all_pos[t_c_end] - 1);
                            new_contog2 = new_fasta_seq.substr(0, all_pos[t_c_start]);

                        } else // transchimera in the middle
                        {
                            new_contig = new_fasta_seq.substr(all_pos[t_c_start], all_pos[t_c_end] - all_pos[t_c_start]);
                            new_contog2 = new_fasta_seq.substr(0, all_pos[t_c_start]) + new_fasta_seq.substr(all_pos[t_c_end], new_fasta_seq.size());
                        }
                        //   cout << "end splitting multigene chimera" << endl;

                        final_pos_file << ID1 << " " << "m_c" << " " << all_pos[t_c_start] << " " << all_pos[t_c_end] << endl;
                        split_contigs.push_back(ID1); // 27oct21
                        // final_pos_file << new_fasta_seq << endl;
                        // final_pos_file << new_contig << endl;


                        //*********************************************************************************************************************                    


                        if (update_ind != string::npos) { //[] found
                            //        cout << "start updating header name" << endl;
                            //string update_str = ID1.substr(update_ind + tool_name_flag.length() + 1, update_ind2 - (update_ind + tool_name_flag.length())); // get everything within :ROAST_ _ROAST,
                            //vector<string> updates;

                            size_t find_BE = ID1.find(ROAST_BE_tag);
                            size_t find_LE = ID1.find(ROAST_LE_tag);
                            size_t find_RE = ID1.find(ROAST_RE_tag);


                            if (find_BE != string::npos) {
                                
                                ID1 = ID1.replace(find_BE, ROAST_LE_tag.size(), ROAST_LE_tag);
                                header1 = ID1.substr(0, update_ind - 1) + "_" + ROAST_newconitg1 + "_" + ID1.substr(update_ind);

                                ID1 = ID1.replace(find_BE, ROAST_RE_tag.size(), ROAST_RE_tag);
                                header2 = ID1.substr(0, update_ind - 1) + "_" + ROAST_newconitg2 + "_" + ID1.substr(update_ind);

                            } else if (find_RE != string::npos) {

                                header1 = ID1.substr(0, update_ind - 1) + "_" + ROAST_newconitg1;
                                header2 = ID1.substr(0, update_ind - 1) + "_" + ROAST_newconitg2 + "_" + tool_name_flag + "_" + ROAST_RE_tag + "_" + tool_name_flag;

                            } else if (find_LE != string::npos) {

                                header1 = ID1.substr(0, update_ind - 1) + "_" + ROAST_newconitg1 + "_" + tool_name_flag + "_" + ROAST_LE_tag + "_" + tool_name_flag;
                                header2 = ID1.substr(0, update_ind - 1) + "_" + ROAST_newconitg2;

                            } else {

                                header1 = ID1.substr(0, update_ind - 1) + "_" + ROAST_newconitg1 + "_" + ID1.substr(update_ind);
                                header2 = ID1.substr(0, update_ind - 1) + "_" + ROAST_newconitg2 + "_" + ID1.substr(update_ind); //ID1.size() - update_ind);

                            }
                        } else { // NO improvement tag yet      

                            header1 = ID1 + "_" + ROAST_newconitg1;
                            header2 = ID1 + "_" + ROAST_newconitg2;
                        }
                        //*****************************************************************************************************************************                  
                        processed_file << ">" << header1 << endl << new_contog2 << endl;
                        processed_file << ">" << header2 << endl << new_contig << endl;
                        split_count++;
                        processed_contigs.push_back(ID_main);
                        ignore_cap3.clear(); // only multi-gene chimeras are fixed not self chimeras so keep their cap3 derivatives

                    } else {

                        if (update_ind != string::npos) { //tag found}

                            //  string update_str = ID1.substr(update_ind + tool_name_flag.length() + 1, update_ind2 - (update_ind + tool_name_flag.length()));
                            size_t find_U = ID1.find(ROAST_update_tag);

                            if (find_U != string::npos) {
                                //ignore
                                header1 = ID1;
                                
                            } else {
                                header1 = ID1.substr(0, ID1.size() - tool_name_flag.length() - 1) + "-" + ROAST_update_tag + "_" + tool_name_flag;
                                update_count++;
                            }
                        } else {
                            header1 = ID1 + "_" + tool_name_flag + "_" + ROAST_update_tag + "_" + tool_name_flag;
                            update_count++;
                        }

                        processed_file << ">" << header1 << endl << new_fasta_seq << endl;
                        processed_contigs.push_back(ID_main);

                    }
                    //  cout << "end updating header name" << endl;
                   /* all_pos.clear();
                     problem_type.clear();
                     same_contig.clear();
                     missing_seq.clear();*/
                }// if(updated)
                all_pos.clear();
                problem_type.clear();
                same_contig.clear();
                missing_seq.clear();
                update = false;

                if (found = temp[1].find("missingSequence") != string::npos)
                    missing_seq.push_back(temp[4]);

                entry1 = entry;
                same_contig.push_back(entry1);

                ID1 = ID2;
                transchimera = false, translocation = false, inv_tra = false, ignore = false, self_chimera = false;
                i = 0;
            }
        }
        temp.clear();
    }


    string ID_temp, entry_fasta, ID;
    ifstream fasta_file;
    fasta_file.open(fastafile.c_str());

    while (!(getline(fasta_file, entry_fasta).eof())) {

        if (entry_fasta[0] == '>') {
            
            utils.str_split(entry_fasta, temp, delimiter);
            ID_temp = entry_fasta;
            ID = entry_fasta.substr(1, entry_fasta.size()); //ID = temp[0].substr(1, temp[0].size() - 1); changed from this at 14/9/2020
            
        } else {//merged_contigs2.find(ID1) == merged_contigs2.end()
            
            if ((find(processed_contigs.begin(), processed_contigs.end(), ID) != processed_contigs.end()) || (find(ignore_cap3.begin(), ignore_cap3.end(), ID) != ignore_cap3.end())) { //current contig is being split or updated
                continue;   
            } else { // when ID didn't find in list of the sequences
                processed_file << ID_temp << endl << entry_fasta << endl;
            }
        }
    }
    final_pos_file.close();
    processed_file.close();
    pot_chimera_file.close();
    fasta_file.close();
    
    cout << "Number of split chimeric contigs in iteration " <<  outer_iteration << ": " << split_count << endl;
    cout << "Number of updated contigs for local mis-assemblies in iteration " <<  outer_iteration << ": " << update_count << endl;
    
    fixed_chimeric_contigs = fixed_chimeric_contigs + split_count;
    updatad_local_misAssemblies  = updatad_local_misAssemblies + update_count;
    return split_count;
}

bool mis_assembly_chimera::compare_positions(vector<int> all_positions, vector<int> problems, int p_start, int p_end, int problem) { // for probelms of start and end position
    utils utils;
    bool is_problem_repeating = false;
    int pos, start, start2, end, end2, insert_pos, insert_pos2;
    vector<int>::iterator it_pos = all_positions.begin();
    
    for (vector<int>::iterator it = problems.begin(); it != problems.end(); it++) {
        if (*it == PROBLEM_MISSING_SEQ) {
            pos = *it_pos;
            it_pos++;
            if (pos >= p_start && pos <= p_end) {
                is_problem_repeating = true;
                break;
            }
        } else if ((*it == PROBLEM_UNSUPPORTED_INSERTION) || (*it == PROBLEM_MIS_ASSEMBLY) || (*it == PROBLEM_SELF_CHIMERA) || (*it == PROBLEM_MULTIGENE_CHIMERA) || (*it == PROBLEM_INVERSION)) {
            start = *it_pos;
            it_pos++;
            end = *it_pos;
            it_pos++;
            
            if(problem == *it && p_start == start && p_end == end) // ignore entry itself
                continue;
            else { //if (problem != *it && (p_start != start || p_end != end)) {
                // start-p_start-p_end-end  or p_start-start-end-p_end  or  p_start-start-p_end-end  or  start-p_start-end-p_end 
                if ((start <= p_start && end >= p_end) || (start >= p_start && end <= p_end) || (start >= p_start && end >= p_end && start < p_end) || (start <= p_start && end <= p_end && p_start < end)) {
                    is_problem_repeating = true;
                    break;
                }
            }
        } else if ((*it == PROBLEM_TRANSLOCATION) || (*it == PROBLEM_INVERSION_TRANSLOCATION)) {
            start = *it_pos;
            it_pos++;
            end = *it_pos;
            it_pos++;
            insert_pos = *it_pos;
            it_pos++;
            
            if(problem == *it && ((p_start == start && p_end == end) || (p_start == insert_pos || p_end == insert_pos)))
                continue;
            else { // if (problem != *it && (p_start != start || p_end != end) && p_start != insert_pos && p_end != insert_pos) { // p_start or p_end should not same as insert positions of other problem  
                // start-p_start-p_end-end  or p_start-start-end-p_end  or  p_start-start-p_end-end  or  start-p_start-end-p_end 
                if ((start <= p_start && end >= p_end) || (start >= p_start && end <= p_end) || (start >= p_start && end >= p_end && start < p_end) || (start <= p_start && end <= p_end && p_start < end)) {
                    is_problem_repeating = true;
                    break;
                }
            }
        }
        else if ((*it == PROBLEM_INVERSION_TRANSLOCATION_TRANSLOCATION) || (*it == PROBLEM_INVERSION_TRANSLOCATION_INVERSION_TRANSLOCATION)) { // six values
            start = *it_pos;
            it_pos++;
            end = *it_pos;
            it_pos++;
            insert_pos = *it_pos;
            it_pos++;
            
            if(problem == *it && ((p_start == start && p_end == end) || (p_start == insert_pos || p_end == insert_pos)))
                continue;
            
            else{ // if (problem != *it && (p_start != start || p_end != end) && p_start != insert_pos && p_end != insert_pos) { // p_start or p_end should not same as insert positions of other problem  
                // start-p_start-p_end-end  or p_start-start-end-p_end  or  p_start-start-p_end-end  or  start-p_start-end-p_end 
                if ((start <= p_start && end >= p_end) || (start >= p_start && end <= p_end) || (start >= p_start && end >= p_end && start < p_end) || (start <= p_start && end <= p_end && p_start < end)) {
                    is_problem_repeating = true;
                    break;
                }
            }
            start2 = *it_pos;
            it_pos++;
            end2 = *it_pos;
            it_pos++;
            insert_pos2 = *it_pos;
            it_pos++;

            if(problem == *it && ((p_start == start2 && p_end == end2) || (p_start == insert_pos2 || p_end == insert_pos2)))
                continue;
            
            else{ //if (problem != *it && (p_start != start2 || p_end != end2) && p_start != insert_pos && p_end != insert_pos2) { // p_start or p_end should not same as insert positions of other problem  
                // start-p_start-p_end-end  or p_start-start-end-p_end  or  p_start-start-p_end-end  or  start-p_start-end-p_end 
                if ((start2 <= p_start && end2 >= p_end) || (start2 >= p_start && end2 <= p_end) || (start2 >= p_start && end2 >= p_end && start2 < p_end) || (start2 <= p_start && end2 <= p_end && p_start < end2)) {
                    is_problem_repeating = true;
                    break;
                }
            }
        } else if ((*it == PROBLEM_INVERSION_TRANSLOCATION_SELFCHIMERA) || (*it == PROBLEM_INVERSION_TRANSLOCATION_MULTIGENECHIMERA)) { // 3 because start2 s same as end and end2 is same insert_pos
            start = *it_pos;
            it_pos++;
            end = *it_pos;
            it_pos++;
            insert_pos = *it_pos;
            it_pos++;
            // start2 = *it_pos;
            it_pos++;
            // end2 = *it_pos; 
            it_pos++;
            
            if(problem == *it && ((p_start == start && p_end == end) || (p_start == insert_pos || p_end == insert_pos)))
                continue;
            
            else { //if (problem != *it && (p_start != start || p_end != end) && p_start != insert_pos && p_end != insert_pos) { // p_start or p_end should not same as insert positions of other problem  
                // start-p_start-p_end-end  or p_start-start-end-p_end  or  p_start-start-p_end-end  or  start-p_start-end-p_end 
                if ((start <= p_start && end >= p_end) || (start >= p_start && end <= p_end) || (start >= p_start && end >= p_end && start < p_end) || (start <= p_start && end <= p_end && p_start < end)) {
                    is_problem_repeating = true;
                    break;
                }
            }
        } else if (*it == PROBLEM_INVERSION_TRANSLOCATION_SELFCHIMERA) { //because insert_pos is same as start2 ans insert_pos2 is same as end
            start = *it_pos;
            it_pos++;
            end = *it_pos;
            it_pos++;
            // insert_pos = *it_pos; 
            it_pos++;
            
            if(problem == *it && p_start == start && p_end == end)
                continue;
            else{ // if (problem != *it && (p_start != start || p_end != end)) { // p_start or p_end should not same as insert positions of other problem  
                // start-p_start-p_end-end  or p_start-start-end-p_end  or  p_start-start-p_end-end  or  start-p_start-end-p_end 
                if ((start <= p_start && end >= p_end) || (start >= p_start && end <= p_end) || (start >= p_start && end >= p_end && start < p_end) || (start <= p_start && end <= p_end && p_start < end)) {
                    is_problem_repeating = true;
                    break;
                }
            }
            start2 = *it_pos;
            it_pos++;
            end2 = *it_pos;
            it_pos++;
            
            if(problem == *it && p_start == start2 && p_end == end2)
                continue;
            else { // if (problem != *it && (p_start != start2 || p_end != end2)) { // p_start or p_end should not same as insert positions of other problem  
                // start-p_start-p_end-end  or p_start-start-end-p_end  or  p_start-start-p_end-end  or  start-p_start-end-p_end 
                if ((start2 <= p_start && end2 >= p_end) || (start2 >= p_start && end2 <= p_end) || (start2 >= p_start && end2 >= p_end && start2 < p_end) || (start2 <= p_start && end2 <= p_end && p_start < end2)) {
                    is_problem_repeating = true;
                    break;
                }
            }
            it_pos++;
        }

    }
    return is_problem_repeating;

}

bool mis_assembly_chimera::compare_positions(vector<int> all_positions, vector<int> problems, int p_start, int p_end, int insert_at, int problem) { // for problems having 3 positions start, end and insert_site
    utils utils;
    bool is_problem_repeating = false;
    int pos, start, start2, end, end2, insert_pos, insert_pos2;
    vector<int>::iterator it_pos = all_positions.begin();
    
    for (vector<int>::iterator it = problems.begin(); it != problems.end(); it++) {
        if (*it == PROBLEM_MISSING_SEQ ) {
            pos = *it_pos;
            it_pos++;
        } else if ((*it == PROBLEM_UNSUPPORTED_INSERTION) || (*it == PROBLEM_MIS_ASSEMBLY) || (*it == PROBLEM_SELF_CHIMERA) || (*it == PROBLEM_MULTIGENE_CHIMERA) || (*it == PROBLEM_INVERSION)) {
            start = *it_pos;
            it_pos++;
            end = *it_pos;
            it_pos++;
            
            if(problem == *it && p_start == start && p_end == end)
                continue;
            
            else{ //if (problem != *it && p_start != start && p_end != end) {
                // start-p_start-p_end(insert_at)-end  or p_start-start-end-p_end  or  p_start-start-p_end-end  or  start-p_start-end-p_end 
                if ((start <= p_start && (end >= p_end || end >= insert_at)) || (start >= p_start && end <= p_end) || (start >= p_start && end >= p_end && start < p_end) || (start <= p_start && end <= p_end && p_start < end)) {
                    is_problem_repeating = true;
                    break;
                }
                if (insert_at > p_end) {
                    if ((start >= p_end && end >= insert_at && start < insert_at)) { //p_end-start-insert_at-end
                        is_problem_repeating = true;
                        break;
                    }
                } else {
                    if (start <= insert_at && end <= p_start && insert_at < end) { // start-insert_at-end-p_start
                        is_problem_repeating = true;
                        break;
                    }
                }
            }
        } else if ((*it == PROBLEM_TRANSLOCATION) || (*it == PROBLEM_INVERSION_TRANSLOCATION)) {
            start = *it_pos;
            it_pos++;
            end = *it_pos;
            it_pos++;
            insert_pos = *it_pos;
            it_pos++;
            
             if(problem == *it && ((p_start == start && p_end == end) || (p_start == insert_pos || p_end == insert_pos)))
                continue;
            
            else{ //if (problem != *it && (p_start != start || p_end != end) && p_start != insert_pos && p_end != insert_pos) { // p_start or p_end should not same as insert positions of other problem  
                if (insert_at > p_end) {
                    // start-p_start-p_end-end  or p_start-start-end-p_end  or  p_start-start-p_end-end  or  start-p_start-end-p_end 
                    if ((start <= p_start && end >= insert_at) || (start >= p_start && end <= insert_at) || (start >= p_start && end >= insert_at && start < insert_at) || (start <= p_start && end <= insert_at && p_start < end)) {
                        is_problem_repeating = true;
                        break;
                    }
                } else {
                    if ((start <= insert_at && end >= p_end) || (start >= insert_at && end <= p_end) || (start >= insert_at && end >= p_end && start < p_end) || (start <= insert_at && end <= p_end && insert_at < end)) {
                        is_problem_repeating = true;
                        break;
                    }
                }

            }
        } else if ((*it == PROBLEM_INVERSION_TRANSLOCATION_TRANSLOCATION) || (*it == PROBLEM_INVERSION_TRANSLOCATION_INVERSION_TRANSLOCATION)) { // six values
            start = *it_pos;
            it_pos++;
            end = *it_pos;
            it_pos++;
            insert_pos = *it_pos;
            it_pos++;
            
             if(problem == *it && ((p_start == start && p_end == end) || (p_start == insert_pos || p_end == insert_pos)))
                continue;
            
            else { // if (problem != *it && (p_start != start || p_end != end) && p_start != insert_pos && p_end != insert_pos) { // p_start or p_end should not same as insert positions of other problem  
                if (insert_at > p_end) {
                    // start-p_start-p_end-end  or p_start-start-end-p_end  or  p_start-start-p_end-end  or  start-p_start-end-p_end 
                    if ((start <= p_start && end >= insert_at) || (start >= p_start && end <= insert_at) || (start >= p_start && end >= insert_at && start < insert_at) || (start <= p_start && end <= insert_at && p_start < end)) {
                        is_problem_repeating = true;
                        break;
                    }
                } else {
                    if ((start <= insert_at && end >= p_end) || (start >= insert_at && end <= p_end) || (start >= insert_at && end >= p_end && start < p_end) || (start <= insert_at && end <= p_end && insert_at < end)) {
                        is_problem_repeating = true;
                        break;
                    }
                }

            }
            start2 = *it_pos;
            it_pos++;
            end2 = *it_pos;
            it_pos++;
            insert_pos2 = *it_pos;
            it_pos++;

             if(problem == *it && ((p_start == start2 && p_end == end2) || (p_start == insert_pos2 || p_end == insert_pos2)))
                continue;
            
            else{ //if (problem != *it && (p_start != start2 || p_end != end2) && p_start != insert_pos2 && p_end != insert_pos2) { // p_start or p_end should not same as insert positions of other problem  
                if (insert_at > p_end) {
                    // start-p_start-p_end-end  or p_start-start-end-p_end  or  p_start-start-p_end-end  or  start-p_start-end-p_end 
                    if ((start2 <= p_start && end2 >= insert_at) || (start2 >= p_start && end2 <= insert_at) || (start2 >= p_start && end2 >= insert_at && start2 < insert_at) || (start2 <= p_start && end2 <= insert_at && p_start < end2)) {
                        is_problem_repeating = true;
                        break;
                    }
                } else {
                    if ((start2 <= insert_at && end2 >= p_end) || (start2 >= insert_at && end2 <= p_end) || (start2 >= insert_at && end2 >= p_end && start2 < p_end) || (start2 <= insert_at && end2 <= p_end && insert_at < end2)) {
                        is_problem_repeating = true;
                        break;
                    }
                }

            }
        } else if ((*it == PROBLEM_INVERSION_TRANSLOCATION_SELFCHIMERA) || (*it == PROBLEM_INVERSION_TRANSLOCATION_MULTIGENECHIMERA)) { // 3 because start2 s same as end and end2 is same insert_pos
            start = *it_pos;
            it_pos++;
            end = *it_pos;
            it_pos++;
            insert_pos = *it_pos;
            it_pos++;
            // start2 = *it_pos;
            it_pos++;
            // end2 = *it_pos; 
            it_pos++;
            
             if(problem == *it && ((p_start == start && p_end == end) || (p_start == insert_pos || p_end == insert_pos)))
                continue;
            
            else { //if (problem != *it && (p_start != start || p_end != end) && p_start != insert_pos && p_end != insert_pos) { // p_start or p_end should not same as insert positions of other problem  
                if (insert_at > p_end) {
                    // start-p_start-p_end-end  or p_start-start-end-p_end  or  p_start-start-p_end-end  or  start-p_start-end-p_end 
                    if ((start <= p_start && end >= insert_at) || (start >= p_start && end <= insert_at) || (start >= p_start && end >= insert_at && start < insert_at) || (start <= p_start && end <= insert_at && p_start < end)) {
                        is_problem_repeating = true;
                        break;
                    }
                } else {
                    if ((start <= insert_at && end >= p_end) || (start >= insert_at && end <= p_end) || (start >= insert_at && end >= p_end && start < p_end) || (start <= insert_at && end <= p_end && insert_at < end)) {
                        is_problem_repeating = true;
                        break;
                    }
                }

            }
        } else if (*it == PROBLEM_INVERSION_TRANSLOCATIONandSELFCHIMERA) { //because insert_pos is same as start2 ans insert_pos2 is same as end
            start = *it_pos;
            it_pos++;
            end = *it_pos;
            it_pos++;
            // insert_pos = *it_pos; 
            it_pos++;
            
             if(problem == *it && p_start == start && p_end == end)
                continue;
            
            else { //if (problem != *it && (p_start != start || p_end != end) && p_start != insert_pos && p_end != insert_pos) { // p_start or p_end should not same as insert positions of other problem  
                if (insert_at > p_end) {
                    // start-p_start-p_end-end  or p_start-start-end-p_end  or  p_start-start-p_end-end  or  start-p_start-end-p_end 
                    if ((start <= p_start && end >= insert_at) || (start >= p_start && end <= insert_at) || (start >= p_start && end >= insert_at && start < insert_at) || (start <= p_start && end <= insert_at && p_start < end)) {
                        is_problem_repeating = true;
                        break;
                    }
                } else {
                    if ((start <= insert_at && end >= p_end) || (start >= insert_at && end <= p_end) || (start >= insert_at && end >= p_end && start < p_end) || (start <= insert_at && end <= p_end && insert_at < end)) {
                        is_problem_repeating = true;
                        break;
                    }
                }

            }
            start2 = *it_pos;
            it_pos++;
            end2 = *it_pos;
            it_pos++;
            
            if(problem == *it && p_start == start2 && p_end == end2)
                continue;
            
            else { //if (problem != *it && (p_start != start || p_end != end) && p_start != insert_pos && p_end != insert_pos) { // p_start or p_end should not same as insert positions of other problem  
                if (insert_at > p_end) {
                    // start-p_start-p_end-end  or p_start-start-end-p_end  or  p_start-start-p_end-end  or  start-p_start-end-p_end 
                    if ((start2 <= p_start && end2 >= insert_at) || (start2 >= p_start && end2 <= insert_at) || (start2 >= p_start && end2 >= insert_at && start2 < insert_at) || (start2 <= p_start && end2 <= insert_at && p_start < end2)) {
                        is_problem_repeating = true;
                        break;
                    }
                } else {
                    if ((start2 <= insert_at && end2 >= p_end) || (start2 >= insert_at && end2 <= p_end) || (start2 >= insert_at && end2 >= p_end && start2 < p_end) || (start2 <= insert_at && end2 <= p_end && insert_at < end2)) {
                        is_problem_repeating = true;
                        break;
                    }
                }

            }
            //insert_pos2 = *it_pos;
            it_pos++;
        }

    }
    return is_problem_repeating;
}

bool mis_assembly_chimera::compare_positions(vector<int> all_positions, vector<int> problems, int position, int problem) { // for problems having one position only
    utils utils;
    bool is_problem_repeating = false;
    int pos, start, start2, end, end2, insert_pos, insert_pos2;
    vector<int>::iterator it_pos = all_positions.begin();
    for (vector<int>::iterator it = problems.begin(); it != problems.end(); it++) {
        if (*it == PROBLEM_MISSING_SEQ) {
            pos = *it_pos;
            it_pos++;
        } else if ((*it == PROBLEM_UNSUPPORTED_INSERTION) || (*it == PROBLEM_MIS_ASSEMBLY) || (*it == PROBLEM_SELF_CHIMERA)) {
            start = *it_pos;
            it_pos++;
            end = *it_pos;
            it_pos++;
            if (position >= start && position <= end) {
                is_problem_repeating = true;
                break;
            }
        } else if (*it == PROBLEM_SELF_CHIMERA) {
            start = *it_pos;
            it_pos++;
            end = *it_pos;
            it_pos++;
        } else if (*it == PROBLEM_INVERSION) {
            start = *it_pos;
            it_pos++;
            end = *it_pos;
            it_pos++;
        } else if ((*it == PROBLEM_TRANSLOCATION) || (*it == PROBLEM_INVERSION_TRANSLOCATION)) {
            start = *it_pos;
            it_pos++;
            end = *it_pos;
            it_pos++;
            insert_pos = *it_pos;
            it_pos++;
        }
        else if ((*it == PROBLEM_INVERSION_TRANSLOCATION_TRANSLOCATION) || (*it == PROBLEM_INVERSION_TRANSLOCATION_INVERSION_TRANSLOCATION)) { // six values
            start = *it_pos;
            it_pos++;
            end = *it_pos;
            it_pos++;
            insert_pos = *it_pos;
            it_pos++;
            start2 = *it_pos;
            it_pos++;
            end2 = *it_pos;
            it_pos++;
            insert_pos2 = *it_pos;
            it_pos++;
        } else if ((*it == PROBLEM_INVERSION_TRANSLOCATION_SELFCHIMERA)) { // 3 because start2 s same as end and end2 is same insert_pos
            start = *it_pos;
            it_pos++;
            end = *it_pos;
            it_pos++;
            insert_pos = *it_pos;
            it_pos++;
            // start2 = *it_pos;
            it_pos++;
            // end2 = *it_pos; 
            it_pos++;
            if (insert_pos > end) {
                if (position >= end && position <= insert_pos) {
                    is_problem_repeating = true;
                    break;
                }
            } else {
                if (position >= insert_pos && position <= start) {
                    is_problem_repeating = true;
                    break;
                }
            }
        } else if (*it == PROBLEM_INVERSION_TRANSLOCATION_MULTIGENECHIMERA) {
            start = *it_pos;
            it_pos++;
            end = *it_pos;
            it_pos++;
            insert_pos = *it_pos;
            it_pos++;
            // start2 = *it_pos;
            it_pos++;
            // end2 = *it_pos; 
            it_pos++;
        } else if (*it == PROBLEM_INVERSION_TRANSLOCATIONandSELFCHIMERA) { //because insert_pos is same as start2 ans insert_pos2 is same as end
            start = *it_pos;
            it_pos++;
            end = *it_pos;
            it_pos++;
            // insert_pos = *it_pos; 
            it_pos++;
            if (position >= end && position <= insert_pos) {
                is_problem_repeating = true;
                break;
                //insert_pos2 = *it_pos;
                it_pos++;
            }
        }
    }
    return is_problem_repeating;

}

void mis_assembly_chimera::update_positions(vector<int> &all_positions, int start, int end, int update_by, bool to_push) {
    for (vector<int>::iterator a = all_positions.begin(); a != all_positions.end(); a++) {
        if (*a > start && *a < end) { //if (*a >= start && *a <= end) {
            if (to_push) { //push
                *a = *a + update_by;
            } else {//if (flag == "pull") { // pull
                *a = *a - update_by;
            }
        }
    }
    //   return all_positions;
    return;
}

void mis_assembly_chimera::get_problemData(vector <string> same_contig, vector<int> &all_pos, vector<int> &problem_type) {
    utils utils;
    //    size_t found;
    //all_pos.clear();
    //problem_type.clear();
    int pos, start_pos, start_pos2, start_pos3, end_pos, end_pos2, end_pos3, insert_at, insert_at2, insert_at3;
    //  multimap <string, vector<int> >all_problems;
    vector<string> temp2, temp;
    int itr = 0;
    while (itr < same_contig.size()) {
        string entry = same_contig[itr];
        utils.str_split(entry, temp, delimiter);
        string problem = temp[1];
        string pos_temp = temp[2];
        if (size_t found = problem.find("missingSequence") != string::npos) {
            pos = atoi(temp[2].c_str());
            all_pos.push_back(pos);
            problem_type.push_back(PROBLEM_MISSING_SEQ);
            //     all_problems.insert(std::pair<string, vector<int> > ("m_s", all_pos)); //insert(make_pair());["m_s"].insert(pos);
            //all_pos.clear();
            //missing_seq.push_back(temp[4]);

        } else if (problem == "miss-assembly") {
            utils.str_split(pos_temp, temp2, "-");
            start_pos = atoi(temp2[0].c_str());
            end_pos = atoi(temp2[1].c_str());
            all_pos.push_back(start_pos);
            all_pos.push_back(end_pos);
            problem_type.push_back(PROBLEM_MIS_ASSEMBLY);
            //     all_problems.insert(std::pair<string, vector<int> > ("m_a", all_pos));
            //all_pos.clear();

        } else if (problem == "unsupported_insertion") {
            utils.str_split(pos_temp, temp2, "-");
            start_pos = atoi(temp2[0].c_str());
            end_pos = atoi(temp2[1].c_str());
            all_pos.push_back(start_pos);
            all_pos.push_back(end_pos);
            problem_type.push_back(PROBLEM_UNSUPPORTED_INSERTION);
            //    all_problems.insert(std::pair<string, vector<int> > ("u_i", all_pos));
            // all_pos.clear();

        } else if (problem == "selfChimera") {
            utils.str_split(pos_temp, temp2, "-");
            start_pos = atoi(temp2[0].c_str());
            end_pos = atoi(temp2[1].c_str());
            all_pos.push_back(start_pos);
            all_pos.push_back(end_pos);
            problem_type.push_back(PROBLEM_SELF_CHIMERA);
            //     all_problems.insert(std::pair<string, vector<int> > ("s_c", all_pos));
            //  all_pos.clear();

        } else if (problem == "multigeneChimera") {
            utils.str_split(pos_temp, temp2, "-");
            start_pos = atoi(temp2[0].c_str());
            end_pos = atoi(temp2[1].c_str());
            all_pos.push_back(start_pos);
            all_pos.push_back(end_pos);
            problem_type.push_back(PROBLEM_MULTIGENE_CHIMERA);
            //    all_problems.insert(std::pair<string, vector<int> > ("t_c", all_pos));
            // all_pos.clear();

        } else if (problem == "translocation") { /// one sided or two sided?
            utils.str_split(pos_temp, temp2, "-");
            start_pos = atoi(temp2[0].c_str());
            end_pos = atoi(temp2[1].c_str());
            string insert_temp = temp[4];
            temp2 .clear();
            utils.str_split(insert_temp, temp2, ":");
            insert_at = atoi(temp2[1].c_str());
            all_pos.push_back(start_pos);
            all_pos.push_back(end_pos);
            all_pos.push_back(insert_at);
            problem_type.push_back(PROBLEM_TRANSLOCATION);
            //     all_problems.insert(std::pair<string, vector<int> > ("t", all_pos));
            //  all_pos.clear();
        } else if (problem == "inversion") { /// one sided or two sided?
            utils.str_split(pos_temp, temp2, "-");
            start_pos = atoi(temp2[0].c_str());
            end_pos = atoi(temp2[1].c_str());
            all_pos.push_back(start_pos);
            all_pos.push_back(end_pos);
            problem_type.push_back(PROBLEM_INVERSION);
            //    all_problems.insert(std::pair<string, vector<int> > ("i", all_pos));
            // all_pos.clear();

        } else if (problem == "inversion_translocation") {
            utils.str_split(pos_temp, temp2, "-");
            start_pos = atoi(temp2[0].c_str());
            end_pos = atoi(temp2[1].c_str());
            string insert_temp = temp[4];
            temp2.clear();
            utils.str_split(insert_temp, temp2, ":");
            insert_at = atoi(temp2[1].c_str());
            all_pos.push_back(start_pos);
            all_pos.push_back(end_pos);
            all_pos.push_back(insert_at);
            problem_type.push_back(PROBLEM_INVERSION_TRANSLOCATION);
            //    all_problems.insert(std::pair<string, vector<int> > ("t_i", all_pos));
            //all_pos.clear();
        } else if (problem == "inversion_translocation_translocation") {
            utils.str_split(pos_temp, temp2, "-");
            start_pos = atoi(temp2[0].c_str());
            end_pos = atoi(temp2[1].c_str());
            string insert_temp = temp[3];
            temp2.clear();
            utils.str_split(insert_temp, temp2, ":");
            insert_at = atoi(temp2[1].c_str());
            pos_temp = temp[4];
            temp2.clear();
            utils.str_split(pos_temp, temp2, "-");
            start_pos2 = atoi(temp2[0].c_str());
            end_pos2 = atoi(temp2[1].c_str());
            insert_temp = temp[5];
            temp2.clear();
            utils.str_split(insert_temp, temp2, ":");
            insert_at2 = atoi(temp2[1].c_str());
            all_pos.push_back(start_pos);
            all_pos.push_back(end_pos);
            all_pos.push_back(insert_at);
            all_pos.push_back(start_pos2);
            all_pos.push_back(end_pos2);
            all_pos.push_back(insert_at2);
            problem_type.push_back(PROBLEM_INVERSION_TRANSLOCATION_TRANSLOCATION);
            //      all_problems.insert(std::pair<string, vector<int> > ("i_t_t", all_pos));
            // all_pos.clear();

        } else if (problem == "inversion_translocation_selfChimera") {
            utils.str_split(pos_temp, temp2, "-");
            start_pos = atoi(temp2[0].c_str());
            end_pos = atoi(temp2[1].c_str());
            string insert_temp = temp[3];
            temp2.clear();
            utils.str_split(insert_temp, temp2, ":");
            insert_at = atoi(temp2[1].c_str());
            pos_temp = temp[4];
            temp2.clear();
            utils.str_split(pos_temp, temp2, "-");
            start_pos2 = atoi(temp2[0].c_str());
            end_pos2 = atoi(temp2[1].c_str());
            all_pos.push_back(start_pos);
            all_pos.push_back(end_pos);
            all_pos.push_back(insert_at);
            all_pos.push_back(start_pos2);
            all_pos.push_back(end_pos2);
            problem_type.push_back(PROBLEM_INVERSION_TRANSLOCATION_SELFCHIMERA);
            //    all_problems.insert(std::pair<string, vector<int> > ("i_t_sc", all_pos));
            // all_pos.clear();

        } else if (problem == "inversion_translocation_inversion_translocation") {
            utils.str_split(pos_temp, temp2, "-");
            start_pos = atoi(temp2[0].c_str());
            end_pos = atoi(temp2[1].c_str());
            string insert_temp = temp[3];
            temp2.clear();
            utils.str_split(insert_temp, temp2, ":");
            insert_at = atoi(temp2[1].c_str());
            pos_temp = temp[4];
            temp2.clear();
            utils.str_split(pos_temp, temp2, "-");
            start_pos2 = atoi(temp2[0].c_str());
            end_pos2 = atoi(temp2[1].c_str());
            insert_temp = temp[5];
            temp2.clear();
            utils.str_split(insert_temp, temp2, ":");
            insert_at2 = atoi(temp2[1].c_str());
            all_pos.push_back(start_pos);
            all_pos.push_back(end_pos);
            all_pos.push_back(insert_at);
            all_pos.push_back(start_pos2);
            all_pos.push_back(end_pos2);
            all_pos.push_back(insert_at2);
            problem_type.push_back(PROBLEM_INVERSION_TRANSLOCATION_INVERSION_TRANSLOCATION);
            //      all_problems.insert(std::pair<string, vector<int> > ("i_t_i_t", all_pos));
            //  all_pos.clear();

        } else if (problem == "inversion_translocation_multigeneChimera") {
            utils.str_split(pos_temp, temp2, "-");
            start_pos = atoi(temp2[0].c_str());
            end_pos = atoi(temp2[1].c_str());
            string insert_temp = temp[3];
            temp2.clear();
            utils.str_split(insert_temp, temp2, ":");
            insert_at = atoi(temp2[1].c_str());
            pos_temp = temp[4];
            temp2.clear();
            //temp2 = utils.split(pos_temp, "-"); // atoi(temp[1].c_str());
            utils.str_split(pos_temp, temp2, "-");
            start_pos2 = atoi(temp2[0].c_str());
            end_pos2 = atoi(temp2[1].c_str());

            all_pos.push_back(start_pos);
            all_pos.push_back(end_pos);
            all_pos.push_back(insert_at);
            all_pos.push_back(start_pos2);
            all_pos.push_back(end_pos2);
            problem_type.push_back(PROBLEM_INVERSION_TRANSLOCATION_MULTIGENECHIMERA);
            //      all_problems.insert(std::pair<string, vector<int> > ("i_t_tc", all_pos));
            // all_pos.clear();
        } else if (problem == "inversion_translocation&selfChimera") {
            utils.str_split(pos_temp, temp2, "-");
            start_pos = atoi(temp2[0].c_str());
            end_pos = atoi(temp2[1].c_str());
            string insert_temp = temp[3];
            temp2.clear();
            utils.str_split(insert_temp, temp2, ":");
            insert_at = atoi(temp2[1].c_str());
            pos_temp = temp[4];
            temp2.clear();
            utils.str_split(pos_temp, temp2, "-"); // atoi(temp[1].c_str());
            start_pos2 = atoi(temp2[0].c_str());
            end_pos2 = atoi(temp2[1].c_str());
            insert_temp = temp[5];
            temp2.clear();
            utils.str_split(insert_temp, temp2, ":");
            insert_at2 = atoi(temp2[1].c_str());
            pos_temp = temp[6];
            temp2.clear();
            utils.str_split(pos_temp, temp2, "-"); // atoi(temp[1].c_str());
            start_pos3 = atoi(temp2[0].c_str());
            end_pos3 = atoi(temp2[1].c_str());

            all_pos.push_back(start_pos);
            all_pos.push_back(end_pos);
            all_pos.push_back(insert_at);
            all_pos.push_back(start_pos2);
            all_pos.push_back(end_pos2);
            all_pos.push_back(insert_at2);
            all_pos.push_back(start_pos3);
            all_pos.push_back(end_pos3);
            problem_type.push_back(PROBLEM_INVERSION_TRANSLOCATIONandSELFCHIMERA);
        }
        itr++;
        temp.clear();
        temp2.clear();
    }
    return;
}
