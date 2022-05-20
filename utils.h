/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   utils.h
 * Author: madiha
 *
 * Created on April 26, 2019, 7:44 PM
 */


/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   processing.h
 * Author: madiha
 *
 * Created on November 5, 2019, 1:12 PM
 */
#include <sstream>
#include <iostream>
#include <fstream>                                                                                                                                                                                         
#include <algorithm>
#include <map>
#include <vector>
#include <list>
#include <set>
#include "filterSamToFastq.h"
//#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <cctype>    
#include <climits>
#include <math.h>       /* ceil */
#include "api/BamAlignment.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "global.h"
#include <boost/algorithm/string/replace.hpp>

using namespace std;
using namespace BamTools;
using namespace boost;

class utils {
public:
   //  int mapping_quality_TH = 20;
    string base;
    string alt;
    std::vector<std::string> indel;
    std::vector<std::string> rID, mID, dir;
    string Rcontig;
    std::vector<std::string> temp, temp2, t;
    string out;
    string lines, linef, ss, s1, s2;

    std::vector<std::string> split(std::string str, std::string sep) {
        char* cstr = const_cast<char*> (str.c_str());
        std::vector<std::string> arr;
        char* token = strtok(const_cast<char*> (str.c_str()), sep.c_str());
        while (token) {
            arr.push_back(token);
            token = strtok(NULL, sep.c_str());
            //token = NULL;
        }
        return arr;
    }

    string Rcomplement(string n) {
        stringstream stream;// = "";
        string str;
        string random = "";
        string reverse_str, header;

        for (int a = 0; a < n.size(); a++) {
            random = n[a];
            if (random.find(">") == 0) {
                header = random;
            } else {
                stream << random;
            }
            random = "";
        }
        str = stream.str();
        stream.str(string()); //clear 
        
        for (int i = 0; i <= str.length(); i++) {
            switch (str[i]) {
                case 'A':
                    str[i] = 'T';
                    break;
                case 'T':
                    str[i] = 'A';
                    break;
                case 'G':
                    str[i] = 'C';
                    break;
                case 'C':
                    str[i] = 'G';
                    break;
            }
        }
        reverse(str.begin(), str.end());
        if (header != "") {
            reverse_str = header + "\n" + str;
            return reverse_str;
        } else return str;
    }

    std::vector<std::string> readfile(string str) {
        ifstream file(str.c_str());
        while (getline(file, lines)) {
            //  if (line[0] == '#') { // ignore header
            //            continue;
            //        }
            temp2.push_back(lines);
        }
        file.close();
        return temp2;
    }

    string fix_SNP(string vcf, string s) {
        std::vector<std::string>snp = readfile(vcf);
        for (int j = 0; j < snp.size(); j++) {
            temp = split(snp[j], "\t");
            base.append(temp[3]);
            alt.append(temp[4]);
            t = split(temp[7], ";");
            indel.push_back(t[0]);
        }
        int count = 0;
        for (int l = 0; l <= indel.size(); l++) {
            if (indel[l] == "INDEL") //ignore indels as they are not present in ref seq
            {
                alt.erase(l, 1);
            }
        }
        for (int n = 0; n < s.length(); n++) {
            if (alt[n] != '.') {
                s[n] = alt[n];
                count = count + 1;
            }
        }
        cout << "No. of SNPs = " << count << endl;
        return s;
    }
    //    string extract_fasta(string ID, string seq_fasta) {
    //        ifstream fastafile(seq_fasta.c_str());
    //        string seq, seq1;
    //        string pre;
    //        string id; //= " ";
    //        int f = 0;
    //        string st;
    //        while (!(getline(fastafile, linef).eof())) {
    //            f = f + 1;
    //            if (linef.find(">") == 0) {
    //                st.clear();
    //                vector<string> temp = split(linef, " ");
    //                string d = temp[0];
    //                string idd = linef.substr(1, linef.size() - 1);
    //                if (idd == ID) {
    //                    st = ID;
    //                    id = linef.substr(1, ID.size());
    //                    seq.clear();
    //                }
    //
    //            } else {
    //                if (id == st) {
    //                    seq += linef;
    //                }
    //            }
    //
    //        }
    //        return seq;
    //    }

    string extract_fasta(string ID, string seq_fasta) {
        ifstream fastafile(seq_fasta.c_str());
        bool ID_found = false;
        string seq, seq1;
        string pre;
        string id; //= " ";
        int f = 0;
        string st;
        while (!(getline(fastafile, linef).eof())) {
            f = f + 1;
            if (linef.find(">") == 0) { //  found
                st.clear();
                vector<string> temp; // = split(linef, " ");
                temp.clear();
                str_split(linef, temp, " ");
                string d = temp[0];
                string idd = d.substr(1, id.size() - 1);
                if (idd == ID) {
                    st = ID;
                    id = linef.substr(1, ID.size());
                    seq.clear();
                    ID_found = true;
                }

            } else { // not found
                if (id == st) {
                    seq += linef;
                }
            }

        }
        fastafile.close();
        if(!ID_found)
        {
            cerr << ID << " not found in " << seq_fasta << endl;
            exit(0);
        }
        return seq;
    }

    void extract_READs_byIDs(string fastq_file_in, string fastq_file_out, string IDs_file) {
        set<string> IDs;
        ifstream fastqfile_in, id_file;
        ofstream fastqfile_out;

        fastqfile_in.open(fastq_file_in.c_str());
        id_file.open(IDs_file.c_str());
        fastqfile_out.open(fastq_file_out.c_str());
        
        bool get_read = false;

        while (!(getline(id_file, linef).eof())) {
            IDs.insert(linef);
        }

        while (!(getline(fastqfile_in, linef).eof())) {
            
            if(linef[0] == '@'){

                string id = linef;
                get_read = false;
                
                if ( IDs.find(id) != IDs.end() ){
                    fastqfile_out << linef << endl; 
                    get_read = true;
                }
                
            }
            else{
                if(get_read){
                    fastqfile_out << linef << endl; 
                }
            }
        }

    }
    
    int split_fasta(string fasta, string out_dir) {

        ofstream sub_fasta_file;
        ifstream fastafile(fasta.c_str());
        string id, all_len, last_entry, pos;
        int length = 0;

        while (!(getline(fastafile, linef).eof())) {

            if (linef.find(">") == 0) {
                vector<string> temp;
                str_split(linef, temp, " ");
                id = temp[0].substr(1, temp[0].size());
                string out_file = out_dir + "/" + id + ".fa";
                sub_fasta_file.open(out_file.c_str());

                sub_fasta_file << ">" << id << endl;
            } else {
                sub_fasta_file << linef << endl;
                sub_fasta_file.close();
            }
        }
        return 0;
    }
    
    int create_gtf(string fasta, string gtf) {
        ofstream gtf_file;
        gtf_file.open(gtf.c_str());
        
        ifstream fastafile(fasta.c_str());
        string id, all_len, last_entry, pos;
        int length = 0;
        
        while (!(getline(fastafile, linef).eof())) {
            if (linef.find(">") == 0) {
                vector<string> temp;
                str_split(linef, temp, " ");
                id = temp[0].substr(1, temp[0].size());
                last_entry = temp[temp.size() -1];
                
                vector<string> temp2;
                str_split(last_entry, temp2, ":");
                all_len = temp2[1];

                vector<string> temp3;
                str_split(all_len, temp3, ",");

                vector<string> temp4;
                str_split(temp3[0], temp4, "-");

                int st = atoi(temp4[0].c_str());
                
                for (int i = 0; i < temp3.size(); i++) {
                    vector<string> temp4;
                    str_split(temp3[i], temp4, "-");
                    int en = atoi(temp4[1].c_str());

                    length = en;
                }
                //gene:ENSG00000000457	superTranscript	exon	1198	3440	0	+	0	transcript_id "gene:ENSG00000000457"	 transcript_id "transcript:ENST00000367770"
                gtf_file << id << "\t" << "superTranscript" << "\t" << "exon" << "\t" << "1" << "\t" << length << "\t" << "." << "\t" << "+" << "\t" << "." << "\t" << "transcript_id " << "\"" << id << "\"" << "\t" << "transcript_id " << "\"" << id << "\"" << endl;
            }
        }
        return 0;
    }

        int create_gff(string fasta, string gff) {
        ofstream gff_file;
        gff_file.open(gff.c_str());
        ifstream fastafile(fasta.c_str());
        string id, all_len, last_entry, pos;
        int length = 0;
        
        while (!(getline(fastafile, linef).eof())) {
            if (linef.find(">") == 0) {
                vector<string> temp;
                str_split(linef, temp, delimiter);
                id = temp[0].substr(1, temp[0].size());
                last_entry = temp[-1];
                
                vector<string> temp2;
                str_split(last_entry, temp2, ":");
                all_len = temp2[1];
                
                vector<string> temp3;
                str_split(all_len, temp3, ",");
                        
                for(int i = 0; i < temp3.size(); i++){
                    vector<string> temp4;
                    str_split(temp[i], temp4, "-");
                    
                    int st = atoi(temp[0].c_str());
                    int en = atoi(temp[1].c_str());
                    
                    length = length + (en - st);
                }
                gff_file << id << "\t" << "." << "\t" << "gene" << "\t" << "1" << "\t" << length << "\t" << "." << "\t" << "." << "\t" << "." << "\t" << "."<< endl;
            }
        }
        return 0;
    }

           int create_gtf_woHeader(string fasta, string gtf) {
        ofstream gtf_file;
        gtf_file.open(gtf.c_str());
        
        ifstream fastafile(fasta.c_str());
        string id, all_len, last_entry, pos;
        int length = 0;
        
        while (!(getline(fastafile, linef).eof())) {
            if (linef.find(">") == 0) {
                vector<string> temp;
                str_split(linef, temp, " ");
                id = linef.substr(1, linef.size());
         
                //gene:ENSG00000000457	superTranscript	exon	1198	3440	0	+	0	transcript_id "gene:ENSG00000000457"	 transcript_id "transcript:ENST00000367770"
                gtf_file << id << "\t" << "superTranscript" << "\t" << "exon" << "\t" << "1" << "\t";
            } else {
                length = linef.size();
                gtf_file << length << "\t" << "." << "\t" << "+" << "\t" << "." << "\t" << "transcript_id " << "\"" << id << "\";" << " " << "gene_id " << "\"" << id << "\";" << endl;
            }
        }
        return 0;
    }
     
    void filter_sam(string bam_in, string bam_out) {
               size_t find_path = bam_in.find(".");
        string temp_bam = bam_in.substr(0, find_path) + "_temp.bam";
        string temp_bam2 = bam_in.substr(0, find_path) + "_temp2.bam";
        string header_sam = bam_in.substr(0, find_path) + "_header.sam";
        
        BamReader reader;
        vector<int> clipsize, read_pos, gen_pos;
        set<string> discard_reads;
        bool pad;

        int count = 0;
        const SamHeader header = reader.GetHeader(); // returns header data object
        // Opens a bam file
        if (!reader.Open(bam_in)) {
            cerr << "Could not open BAM file." << endl;
            exit(0);
        }

        vector<RefData> ref;
        ref = reader.GetReferenceData();
        BamWriter writer;
        if (!writer.Open(temp_bam, header, ref)) {
            cerr << "Could not open output BAM file" << endl;
            exit(0);
            // return;
        }
        BamTools::BamAlignment al;
        while (reader.GetNextAlignment(al)) {
            if (!(al.IsMateMapped()) || !(al.IsMapped())) { //if read unmapped
                writer.SaveAlignment(al);
                            
            } else if (al.MapQuality >= mapping_quality_TH &&  discard_reads.find(al.Name) == discard_reads.end() && al.IsPrimaryAlignment())  { // to avoid secondary and supplementary alignments
                
               // cout << al.Name << " " << al.MapQuality << endl;
                if (al.RefID == al.MateRefID) { //map in same contigs

                    if ((al.IsReverseStrand()) && (al.IsMateReverseStrand())) { // 00 both on reverse strand // -> ->
                        continue;
                    } else if (!(al.IsReverseStrand()) && !(al.IsMateReverseStrand())) { // 11 both on forward strand <- <-
                        continue;
                    }// ignore if reads are mapped in the wrong direction (<- -> instead of -> <-)
                    else if (!(al.IsReverseStrand()) && al.InsertSize < 0) // ((flag & FLAG_READ_REVERSE_STRAND) && read_base_pos <= mate_base_pos) {
                    {
                        continue;
                    } else if ((al.IsReverseStrand()) && al.InsertSize > 0) {
                        continue;
                    } else {
                        writer.SaveAlignment(al);
                    }
                } else {
                    writer.SaveAlignment(al);
                }
            } else{
                discard_reads.insert(al.Name);
                continue;
                }
            
        }        reader.Close();
        writer.Close();

        
      //  cout << "read temp_bam again to remove reads present in discard_reads and write in out_bam" << endl;
        if (!reader.Open(temp_bam)) {
            cerr << "Could not open BAM file." << endl;
            exit(0);
        }
        if (!writer.Open(temp_bam2, header, ref)) {
            cerr << "Could not open output BAM file" << endl;
            exit(0);
            // return;
        }
        
        while (reader.GetNextAlignment(al)) {
            if (discard_reads.find(al.Name) == discard_reads.end()) {

                writer.SaveAlignment(al);

            } else {
                continue;
            }

        }
        reader.Close();
        writer.Close();
        
        string extract_samHeader = "samtools view -H  " + temp_bam2 + " > " + header_sam;
        std::system(extract_samHeader.c_str());
        
        string bam_reheader = "samtools reheader  " +  header_sam + " " + temp_bam2 + " > " + bam_out;
        std::system(bam_reheader.c_str());
        
        remove(temp_bam.c_str());
        remove(temp_bam2.c_str());
        remove(header_sam.c_str());
    }

    void filter_sam_distant(string bam_file, string bam_out) {
        
        size_t find_path = bam_file.find(".");
        string temp_bam = bam_file.substr(0, find_path) + "_temp.bam";
        
        BamReader reader;
        bool pad;
        
        float one_side_SC = float(one_side_allowed_SCs_RI /read_length)*100;
        float each_side_SC = float(each_side_allowed_SCs_RI/read_length)*100;

        set<string> discard_reads;
        
        vector<int> clipsize, read_pos, gen_pos;
        vector<CigarOp> cigar;
        const SamHeader header = reader.GetHeader(); // returns header data object
        // Opens a bam file
        if (!reader.Open(bam_file)) {
            cerr << "Could not open BAM file." << endl;
            exit(0);
        }
        // RefData read_id, mate_id;
        vector<RefData> ref;
        ref = reader.GetReferenceData();
        BamWriter writer;
        if (!writer.Open(temp_bam, header, ref)) {
            cerr << "Could not open output BAM file" << endl;
            exit(0);
            // return;
        }
        BamTools::BamAlignment al;
        while (reader.GetNextAlignment(al)) {

            //  if (al.MapQuality > 0 ) { // don't consider reads of 0 quality

            if (al.RefID != -1 && al.MateRefID != -1 && discard_reads.find(al.Name) == discard_reads.end()) { // to avoid unaligned reads which gives -1 RefID

                RefData read_id = ref.at(al.RefID);
                RefData mate_id = ref.at(al.MateRefID);

                if (!(al.IsMapped())) {

                    continue; //writer.SaveAlignment(al); //write for un mapped reads

                } else if ((mate_id.RefName.compare(read_id.RefName) == 0) || (mate_id.RefName.compare("*") == 0)) { //map on same contig

                    continue;

                } else if (al.Position <= left_edge_boundary && al.IsReverseStrand() && ((al.MatePosition <= left_edge_boundary && al.IsMateReverseStrand()) || (al.MatePosition >= mate_id.RefLength - right_edge_boundary - 1 && (!(al.IsMateReverseStrand()))))) { //in rev direction

                    if (al.GetSoftClips(clipsize, read_pos, gen_pos) == true) { // consider only those reads which don't have Softclips at internal or both sides
                       
                        if (clipsize.size() == 1) { // one side SC
                            
                            if (clipsize[0] != read_pos[0] && clipsize[0] <= each_side_SC) { // right side soft clips > 12% skip reads

                                writer.SaveAlignment(al);

                            } else if (clipsize[0] == read_pos[0] && clipsize[0] <= one_side_SC) { // left side SCs skip if > 25% because they can map in other contigs

                                writer.SaveAlignment(al);

                            } else if (clipsize.size() > 1 && clipsize[0] <= each_side_SC && clipsize[1] <= each_side_SC) { // or both side soft clips > 12% skip reads

                                writer.SaveAlignment(al);

                            } else {

                                discard_reads.insert(al.Name);
                                //continue;
                            }
                        } else { // both side SCs
                            if (clipsize[0] <= each_side_SC && clipsize[1] <= each_side_SC) { // or both side soft clips > 12% skip reads
                                writer.SaveAlignment(al);
                            } else {

                                discard_reads.insert(al.Name);
                                //continue;
                            }
                        }
                        clipsize.clear();
                        read_pos.clear();
                        gen_pos.clear();
                        
                    } else {
                        writer.SaveAlignment(al);
                    }
                } else if (((al.Position >= read_id.RefLength - right_edge_boundary - 1) && !(al.IsReverseStrand())) && ((al.MatePosition <= left_edge_boundary && al.IsMateReverseStrand()) || (al.MatePosition >= mate_id.RefLength - right_edge_boundary - 1 && (!(al.IsMateReverseStrand()))))) { // in for direction  

                    if (al.GetSoftClips(clipsize, read_pos, gen_pos) == true) { // consider only those reads which don't have Softclips at internal or both sides

                        if (clipsize.size() == 1) { // one side SC
                            if (clipsize[0] == read_pos[0] && clipsize[0] <= each_side_SC) { // left side or both side soft clips -> ignore reads

                                writer.SaveAlignment(al);

                            } else if (clipsize[0] != read_pos[0] && clipsize[0] <= one_side_SC) { // left side SCs skip if > 25% because they can map in other contigs

                                writer.SaveAlignment(al);

                            } else {
                                discard_reads.insert(al.Name);
                                //continue;
                            }
                        } else { // both sides SCs
                            if (clipsize[0] < each_side_SC && clipsize[1] < each_side_SC) {
                                writer.SaveAlignment(al);
                            } else {
                                discard_reads.insert(al.Name);
                                //continue;
                            }
                        }
                        clipsize.clear();
                        read_pos.clear();
                        gen_pos.clear();

                    } else {
                        writer.SaveAlignment(al);
                    }
                } else {
                    discard_reads.insert(al.Name);
                   // continue;
                }
                // }
            }// else { // add unaligned reads which gives -1 RefID
            //   writer.SaveAlignment(al);
          //  }


        }

        reader.Close();
        writer.Close();

        //  cout << "read temp_bam again to remove reads present in discard_reads and write in out_bam" << endl;
        if (!reader.Open(temp_bam)) {
            cerr << "Could not open BAM file." << endl;
            exit(0);
        }
        if (!writer.Open(bam_out, header, ref)) {
            cerr << "Could not open output BAM file" << endl;
            exit(0);
            // return;
        }

        while (reader.GetNextAlignment(al)) {
            if (discard_reads.find(al.Name) == discard_reads.end()) {

                writer.SaveAlignment(al);

            } else {
                continue;
            }

        }
        reader.Close();
        writer.Close();
        remove(temp_bam.c_str());

    }

     void filter_sam_unmapped(string bam_file, string bam_out) {
        BamReader reader;
        bool pad;
        //cout << read_length << " " << one_side_allowed_SCs_RI <<  " " << each_side_allowed_SCs_RI << endl;
        float one_side_SC = float (float(one_side_allowed_SCs_RI) /float(read_length) ) * 100;
        float each_side_SC = float(float(each_side_allowed_SCs_RI) /float(read_length) ) * 100;
        
       // cout << one_side_SC << " " << each_side_SC << endl;
        
      //  int FLAG_READ_REVERSE_STRAND = 16;
       // int FLAG_MATE_REVERSE_STRAND = 32;
       // int FLAG_READ_UNMAPPED = 4;
      //  int FLAG_MATE_UNMAPPED = 8;
        
        int count = 0, len_read_contig;
       // unmapped_contig.reserve(contig_count);
        
        string contig_name, prev_contig = "";
        bool contig_written = false;
        
        vector<int> clipsize, read_pos, gen_pos;
        vector<CigarOp> cigar;

        
        const SamHeader header = reader.GetHeader(); // returns header data object
        // Opens a bam file
        if (!reader.Open(bam_file)) {
            cerr << "Could not open BAM file." << endl;
            exit(0);
        }
        // RefData read_id, mate_id;
        vector<RefData> ref;
        ref = reader.GetReferenceData();
        BamWriter writer;
        if (!writer.Open(bam_out, header, ref)) {
            cerr << "Could not open output BAM file" << endl;
            exit(0);
            // return;
        }
        BamTools::BamAlignment al;
        while (reader.GetNextAlignment(al)) {

            //  if (al.MapQuality > 0 ) { // don't consider reads of 0 quality

            if (al.RefID != -1) { // to avoid unaligned reads which gives -1 RefID

                RefData read_id = ref.at(al.RefID);
               // RefData mate_id = ref.at(al.MateRefID);
                len_read_contig = read_id.RefLength;
                contig_name = read_id.RefName;

                
                if (!(al.IsMapped())) {

                        writer.SaveAlignment(al); //write for un mapped reads
                        
                } else if (!(al.IsMateMapped()) && (al.MapQuality > 0)) {

                    if ((al.Position <= left_edge_boundary) && (al.IsReverseStrand())) { //in rev direction

                        if (al.GetSoftClips(clipsize, read_pos, gen_pos) == true) { // consider only those reads which don't have Softclips at internal or both sides

                            if ((clipsize[0] != read_pos[0] && clipsize[0] <= one_side_SC) || (clipsize.size() > 1 && clipsize[0] <= each_side_SC && clipsize[1] < 3)) { // right side or both side soft clips -> ignore reads inner side <3
                                writer.SaveAlignment(al);
                                contig_written = true;

                                unmappded_count ++;

                            } else {
                                continue;
                            }
                            clipsize.clear();
                            read_pos.clear();
                            gen_pos.clear();
                        } else {
                            writer.SaveAlignment(al);
                            contig_written = true;

                            unmappded_count++;
                        }
                    } else if ((al.Position >= read_id.RefLength - right_edge_boundary - 1) && (!(al.IsReverseStrand()))) { // in for direction  

                        if (al.GetSoftClips(clipsize, read_pos, gen_pos) == true) { // consider only those reads which don't have Softclips at internal or both sides
                           
                            if ((clipsize[0] == read_pos[0] && clipsize[0] <= one_side_SC) || (clipsize.size() > 1 && clipsize[0] < 3 && clipsize[1] <= each_side_SC)) { // left side or both side soft clips -> ignore reads, inner side < 3
                                writer.SaveAlignment(al);
                                contig_written = true;

                                unmappded_count++;
                            } else {
                                continue;
                            }

                            clipsize.clear();
                            read_pos.clear();
                            gen_pos.clear();
                        } else {
                            writer.SaveAlignment(al);
                            contig_written = true;

                            unmappded_count++;
                        }
                    }

                }
            } else { // add unaligned reads which gives -1 RefID
                writer.SaveAlignment(al);
            }
            //}
            if(contig_written == true && contig_name != prev_contig){
                unmapped_contig.push_back(contig_name);
                prev_contig = contig_name;
                contig_written = false;
            }
        }

        reader.Close();
        writer.Close();
    }

     
    void filter_sam_extension(string sam_in, string sam_out, string fasta) { //filter sam file to minimize time for next softclip extension
        int FLAG_READ_REVERSE_STRAND = 16;
        int FLAG_MATE_REVERSE_STRAND = 32;
        int FLAG_READ_UNMAPPED = 4;
        int FLAG_MATE_UNMAPPED = 8;
        int length;
        /*save ID and length of contig from fasta file*/
        map<string, int> id_length;
        string linef, ID, str_length = "";

        ifstream fastafile(fasta.c_str());
        while (!(getline(fastafile, linef).eof())) {
            if (linef[0] == '>') {
                //                std::size_t found = linef.find("path");  // for normal fasta file use this
                vector<string> temp = split(linef, " ");
                ID = temp[0].substr(1, temp[0].length() - 1);
                //                string two = temp[1];
                //                length = atoi(two.substr(4, two.length() - 1).c_str());
                //                id_length[ID] = length;
                continue;
            } else {
                length = linef.size(); //for supertranscript use this
                id_length[ID] = length;

            }
        }
        fastafile.close();

        string entry, read_contig, mate_contig, cigar, unmapped_read, unmapped_entry, read_name;
        int flag, read_base_pos, mate_base_pos, insert_size, mapping_quality;

        vector<string> temp, added_read;
        ofstream out_file;
        out_file.open(sam_out.c_str());
        ifstream sam_file(sam_in.c_str());
        while (!(getline(sam_file, entry).eof())) {
            int softclip = 0;
            if (entry[0] == '@') {
                out_file << entry << endl;
            } else {
                temp = split(entry, "\t");
                read_base_pos = atoi(temp[3].c_str());
                read_name = temp[0];
                read_contig = temp[2];
                mate_contig = temp[6];
                mate_base_pos = atoi(temp[7].c_str());
                insert_size = atoi(temp[8].c_str());
                mapping_quality = atoi(temp[4].c_str());
                flag = atoi(temp[1].c_str());
                cigar = temp[5];

                for (int i = 0; i < cigar.length(); i++) {//check 
                    if (isdigit(cigar[i])) {
                        str_length += cigar[i];
                    } else if (cigar[i] == 'S') {
                        softclip = softclip + atoi(str_length.c_str());
                        str_length = "";
                    } else if (cigar[i] == 'D') {
                        str_length = "";
                    } else if (cigar[i] == 'I') {
                        str_length = "";
                    } else if (cigar[i] == 'M') {
                        str_length = "";
                    }
                }
                int len_read_contig = id_length[read_contig];
                // if (len_read_contig < 300 || softclip > 10) continue;
                //                if (mapping_quality < mapping_quality_TH) //|| softclip > 10) 
                //                    continue;
                //   else {'
                if ((read_base_pos <= 100) || (read_base_pos >= len_read_contig - 100) || (find(added_read.begin(), added_read.end(), read_name) < added_read.end())) {
                    added_read.push_back(read_name);
                    out_file << entry << endl;
                }
            }

        }
    }

    void embeded_seq(string blastn_file, string fastafile, string new_assembly, string embeded_contigs) {
        vector <string> embeded_seq, temp;
        string contig1, contig2, entry_blast, entry_fasta, ID;
        float p_identity, aln_len, query_len;
        bool embeded_found = false;
        ifstream blastn, f_file;
        blastn.open(blastn_file.c_str());
        ofstream newAssembly;
        ofstream embeded_assembly;

        while (!(getline(blastn, entry_blast).eof())) {
            if (entry_blast[0] == '#') {
                embeded_found = false;
                continue;
            } else if (embeded_found == true) { // if embeded found once for one contig
                continue;
            } else {
                temp = split(entry_blast, "\t");
                contig1 = temp[0];
                contig2 = temp[1];
                p_identity = atoi(temp[2].c_str());
                aln_len = atoi(temp[12].c_str());
                query_len = atoi(temp[8].c_str());
                if (contig1.substr(0, contig1.size() - 3) == contig2.substr(0, contig2.size() - 3)) { //BOTH CONTIGS FROM SAME CLUSTER AN D HAVE SAME NAME, ignore them
                    embeded_found = false;
                    continue;
                } else if (find(embeded_seq.begin(), embeded_seq.end(), contig2) != embeded_seq.end()) { //contig2 is already added as contig1 in embedded
                    embeded_found = false;
                    continue;
                } else {
                    float query_cov = (aln_len / query_len)*100;
                    if (query_cov >= 90 && p_identity >= 80) { // check identity and query coverage
                        embeded_seq.push_back(contig1);
                        embeded_found = true;
                    } else {
                        embeded_found = false;
                        continue;
                    }
                }
            }

        } //end reading blastn file
        blastn.close();

        f_file.open(fastafile.c_str());
        newAssembly.open(new_assembly.c_str());
        embeded_assembly.open(embeded_contigs.c_str());
        while (!(getline(f_file, entry_fasta).eof())) {
            if (entry_fasta[0] == '>') {
                temp = split(entry_fasta, "\t");
                ID = temp[0].substr(1, temp[0].size() - 1);
            } else {
                if (find(embeded_seq.begin(), embeded_seq.end(), ID) == embeded_seq.end()) {
                    newAssembly << ">" << ID << endl << entry_fasta << endl;
                } else { // when ID found in list of embeded sequnces
                    embeded_assembly << ">" << ID << endl << entry_fasta << endl;
                }


            }


        }
        newAssembly.close();
        embeded_assembly.close();
    }

    void blast_analysis(string blastn_file, string analysis_result_file) { // not whole assembly blast instead blast only identified contigs
        ifstream blastn;
        string entry_blast, contig1, contig2, ID;
        float aln_len, query_len, p_identity;
        bool check = false;
        vector<string> temp;
        float count = 0.0, total_merged = 0.0;

        blastn.open(blastn_file.c_str());
        ofstream out_file;
        out_file.open(analysis_result_file.c_str());
        while (!(getline(blastn, entry_blast).eof())) {
            //temp = split(entry_blast, " ");
            temp.clear();
            str_split(entry_blast, temp, " ");
            if ((temp[0] == "#") && (temp[1] == "Query:")&& (temp[3] == "&")) {
                //temp = split(entry_blast, ":");
                str_split(entry_blast, temp, ":");
                ID = temp[1];
                check = false;
                total_merged++;
                continue;
            } else if (entry_blast[0] == '#') {
                check = false;
                continue;
            } else if (check == true) { // if embeded found once for one contig
                continue;
            } else if (ID != "") {
                //temp = split(entry_blast, "\t");
                str_split(entry_blast, temp, delimiter);
                contig1 = temp[0];
                // contig2 = temp[1]; //for blastx
                contig2 = temp[2];
                //  p_identity = atof(temp[5].c_str());
                // aln_len = atof(temp[7].c_str());
                // query_len = atof(temp[1].c_str());
                //   float query_cov = (aln_len / query_len)*100;
                //  float query_cov = atof(temp[5].c_str()); //for blastx //for drosophilla
                float query_cov = atof(temp[4].c_str());
                // float  q_cov =  ((aln_len) / query_len)*100;
                if (query_cov >= 70) {
                    check = true;
                    out_file << ID << endl << entry_blast << endl;
                    ID = "";
                    count++;

                } else {
                    check = false;
                    continue;
                }
            }

        } //end reading blastn file
        blastn.close();
        out_file.close();
        cout << "number of confirmed merged contigs : " << count << endl;
        cout << "number of total merged contigs : " << total_merged << endl;
        cout << "percentage of verified merged_contigs: " << (count / total_merged)*100 << endl;
    }

    void blast_individual_analysis(string blastn_file, string analysis_result_file, string merged_file) {// change last input file***
        ifstream blastn;
        string entry_blast, contig1, contig2, ID, entry_fasta;
        float aln_len, query_len, p_identity;
        bool check = false;
        vector<string> temp, temp1, temp2;
        float count = 0.0, total_merged = 0.0;
        int id = 0;
        map<string, string> merged_header, blast_result;
        ifstream fastafile;
        fastafile.open(merged_file.c_str());



        typedef multimap<string, string > mmap;
        mmap merged;
        typedef pair<string, string> id1_id2;
        ptrdiff_t itt; // pointer that get match

        /* when check all broken */
        while (!(getline(fastafile, entry_fasta).eof())) { //***
            temp1 = split(entry_fasta, "\t");
            total_merged++;
            merged_header.insert(make_pair(temp1[0], temp1[4])); // read and mate contig
        }
        fastafile.close();

        /* when check merged file*/
        //        while (!(getline(fastafile, entry_fasta).eof())) { //***
        //                temp1 = split(entry_fasta,";");
        //                temp2 = split(temp1[1], ",");
        //                total_merged++;
        //                merged_header.insert(make_pair(temp2[0], temp2[1])); // read and mate contig
        //        }
        //        fastafile.close();

        blastn.open(blastn_file.c_str());
        ofstream out_file;
        out_file.open(analysis_result_file.c_str());
        while (!(getline(blastn, entry_blast).eof())) {
            temp = split(entry_blast, " ");
            if ((temp[0] == "#") && (temp[1] == "Query:")) {
                check = false;
                id = 0;
                continue;
            } else if (entry_blast[0] == '#') {
                id = 0;
                check = false;
                continue;
            } else if (check == true) { // if embeded found once for one contig
                id = 0;
                continue;
            } else if (id < 3) {
                temp = split(entry_blast, "\t");
                contig1 = temp[0];
                contig2 = temp[2];
                //  blast_result.insert(make_pair(contig1, contig2)); // read and mate contig
                merged.insert(id1_id2(contig1, contig2));
                id++;
                check = false;
                continue;
            } else {
                check = true;
                id = 0;
            }

        } //end reading blastn file

        vector<string> subject1;
        string subject2 = "";
        int s2 = 0;

        for (map<string, string>::iterator it_fasta = merged_header.begin(); it_fasta != merged_header.end(); ++it_fasta) {
            std::pair <std::multimap<string, string>::iterator, std::multimap<string, string>::iterator> ret1;
            ret1 = merged.equal_range(it_fasta->first); //extract contig1 Hits from blast data map and get 1 or top 3 hits as user define
            for (std::multimap<string, string>::iterator it = ret1.first; it != ret1.second; ++it) {

                subject1.push_back(it->second); //  put all hits for contig1 in subject1
            }
            std::pair <std::multimap<string, string>::iterator, std::multimap<string, string>::iterator> ret2;
            ret2 = merged.equal_range(it_fasta->second); // extract contig2 Hits from blast data map and get 1 or top 3 hits as user define

            for (std::multimap<string, string>::iterator it = ret2.first; it != ret2.second; ++it) {
                itt = std::find(subject1.begin(), subject1.end(), it->second) - subject1.begin(); // search hit of contig2 in contig1(subject1)
                if (itt != subject1.size()) { // if match found
                    subject2 = it->second;
                    out_file << it_fasta->first << "\t" << subject1[itt] << "\t" << it_fasta->second << "\t" << subject2 << "\t" << itt << "\t" << s2 << endl;
                    count++;
                    break;
                }
                s2++;
            }

            subject1.clear();
            subject2 = "";
            itt = 0;
            s2 = 0;

        }


        blastn.close();
        out_file.close();
        cout << "number of confirmed merged contigs : " << count << endl;
        cout << "number of total merged contigs : " << total_merged << endl;
        cout << "percentage of verified merged_contigs: " << (count / total_merged)*100 << endl;
    }


    //    void broken_checkby_BLAST(string blastn_file, string analysis_result_file, string broken_file) { // to check FP and FN from over all data // for now just from Hits of RI
    //        ifstream blastn;
    //        string entry_blast, contig1="", contig2="", ID, entry_fasta;
    //        float aln_len, query_len, p_identity;
    //        bool check = false;
    //        vector<string> temp, temp1, temp2;
    //        float count = 0.0, total_merged = 0.0;
    //        int id = 0, max_read_no = 0;
    //        map<string, string> merged_header, blast_result;
    //         std::map<string, string>::iterator it;
    //        ifstream fastafile;
    //        fastafile.open(broken_file.c_str());
    //
    //
    //
    //        typedef multimap<string, string > mmap;
    //        mmap blast_merged;
    //        typedef pair<string, string> id1_id2;
    //        ptrdiff_t itt; // pointer that get match
    //
    //
    //        while (!(getline(fastafile, entry_fasta).eof())) {
    //            if(entry_fasta[0]=='r') continue; //ignore first line
    //            else{
    //                temp1 = split(entry_fasta,"\t");
    //                int read_no = atoi(temp1[1].c_str());
    //             //  it = merged_header.find(temp1[4]);
    //
    ////                if( temp1[4] == "TRINITY_DN10506_c0_g1") {
    ////                cout<<temp1[0] << "   " <<temp1[4] <<endl;
    ////            }  
    //                if(temp1[0].substr(0,temp1[0].size()-6) == temp1[4].substr(0,temp1[4].size()-6)){
    //                 continue; // to ignore same gene id names of clusters              
    //                }
    //                else if((temp1[0] == contig1) && ( read_no > max_read_no)){ // to get contig2 with max number of reads for multiple hits of contig1
    //                contig1 = temp1[0];
    //                contig2 = temp1[4];
    //                max_read_no = read_no;
    //                }
    //                
    //                else if ((temp1[0] != contig1) && (merged_header.find(temp1[4]) == merged_header.end())){ // to get unique entr  y
    //                merged_header.insert(make_pair(contig1, contig2)); // read and mate contig
    //                }
    //                else{
    //                    continue;
    //                }
    //                contig1 = temp1[0];
    //                contig2 = temp1[4];
    //                max_read_no = atoi(temp1[1].c_str());
    //                  }
    //        }
    //        fastafile.close();
    //
    //        blastn.open(blastn_file.c_str());
    //        ofstream out_file;
    //        out_file.open(analysis_result_file.c_str());
    //        while (!(getline(blastn, entry_blast).eof())) {
    //            temp = split(entry_blast, " ");
    //            if ((temp[0] == "#") && (temp[1] == "Query:")) {
    //                check = false;
    //                id = 0;
    //                continue;
    //            } 
    //            else if (entry_blast[0] == '#') {
    //                id = 0;
    //                check = false;
    //                continue;
    //            } 
    //            else if (check == true) { // if embeded found once for one contig
    //                id = 0;
    //                continue;
    //            }
    //            
    //            else if (id < 3) {
    //                temp = split(entry_blast, "\t");
    //                contig1 = temp[0];
    //                contig2 = temp[2];
    //              //  blast_result.insert(make_pair(contig1, contig2)); // read and mate contig
    //                blast_merged.insert(id1_id2(contig1, contig2));
    //                id++;
    //                check = false;
    //                continue;
    //            } else {
    //                check = true;
    //                id = 0;
    //            }
    //
    //        } //end reading blastn file
    //
    //        vector<string> subject1;
    //        string subject2 = ""; 
    //        int s2 =0;
    //
    //        for (map<string, string>::iterator it_fasta = merged_header.begin(); it_fasta != merged_header.end(); ++it_fasta) {
    //            total_merged++;
    //            std::pair <std::multimap<string, string>::iterator, std::multimap<string, string>::iterator> ret1;
    //            ret1 = blast_merged.equal_range(it_fasta->first);  //extract contig1 Hits from blast data map and get 1 or top 3 hits as user define
    //            for (std::multimap<string, string>::iterator it = ret1.first; it != ret1.second; ++it) {
    //
    //                subject1.push_back(it->second); //  put all hits for contig1 in subject1
    //            }
    //            std::pair <std::multimap<string, string>::iterator, std::multimap<string, string>::iterator> ret2;
    //            ret2 = blast_merged.equal_range(it_fasta->second); // extract contig2 Hits from blast data map and get 1 or top 3 hits as user define
    //            
    //            for (std::multimap<string, string>::iterator it = ret2.first; it != ret2.second; ++it) {                                            
    //                itt = std::find(subject1.begin(), subject1.end(), it->second) - subject1.begin(); // search hit of contig2 in contig1(subject1)
    //                if (itt != subject1.size()) {  // if match found
    //                    subject2 = it->second;
    //                    out_file << it_fasta->first << "\t" << subject1[itt] << "\t" << it_fasta->second << "\t" << subject2 << "\t" << itt << "\t"<< s2 <<endl;
    //                    count++; 
    //                    break;
    //                }
    //                else
    //                    cout<< it_fasta->first<< "\t" << it_fasta->second<<endl;
    //                s2 ++;
    //            }
    //
    //            subject1.clear();
    //            subject2 = "";
    //            itt = 0; s2 = 0;
    //
    //        }
    //
    //
    //        blastn.close();
    //        out_file.close();
    //        cout << "number of confirmed merged contigs : " << count << endl;
    //        cout << "number of total merged contigs : " << total_merged << endl;
    //        cout << "percentage of verified merged_contigs: " << (count / total_merged)*100 << endl;
    //    }

    void broken_checkby_BLAST(string blastn_file, string analysis_result_file) { // to check errors in over all data 
        string entry;
        bool check = true;
        int query_cov, query_cov_TH = 85, sub_st, sub_en, sub_st1, sub_en1, id = 0;
        float q_len, sub_len;
        float q_st, q_end, q_st1, q_end1;
        
        vector<string> temp;
        map<string, vector<string> > entry_blast;
        vector<string> blast_temp, blast_temp2;
        map<string, vector<string> >::iterator entries_it, it;
        // coverage[temp3[0]].insert (make_pair (pos,cov)); 
        ofstream blast_analysis;
        blast_analysis.open(analysis_result_file.c_str());
        ifstream blastn;
        blastn.open(blastn_file.c_str());

        /*ALl possible Broken contigs*/
        /*# Fields: query id, query length, subject id, subject length, % query coverage per subject, query/sbjct frames,
         *  % identity, identical, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue,
         *  bit score */
        while (!(getline(blastn, entry).eof())) {
            temp = split(entry, "\t");
            if (entry[0] == '#') {
                id = 0;
                check = false;
                continue;
            } else if (check == true) { // if embeded found once for one contig
                id = 0;
                continue;
            } else if (id < 1) {
                temp = split(entry, "\t");
                query_cov = atoi(temp[4].c_str());
                vector<string> temp2;
                vector<string> temp_dup_check;
                temp_dup_check.push_back(temp[2]); // save subject name
                temp_dup_check.push_back(temp[4]); // save query coverage
                if (id == 1) {

                }
                q_len = atof(temp[8].c_str());
                sub_len = atof(temp[3].c_str());
                float q_c = float((q_len / sub_len)*100);
                //if (q_c >= 40) 
                temp2.push_back(temp[0]); // query_id
                temp2.push_back(temp[13]); // s_start
                temp2.push_back(temp[14]); //s_end
                temp2.push_back(temp[11]); //q_start
                temp2.push_back(temp[12]); // q_end
                temp2.push_back(temp[1]); // query_length
                temp2.push_back(temp[3]); // subject_length

                if (query_cov >= 0 && q_c >= 0) {// (query_cov >= query_cov_TH && q_c >= 10)
                    entries_it = entry_blast.find(temp[2]); // another query found for same subject
                    if (entries_it == entry_blast.end())
                        entry_blast.insert(make_pair(temp[2], temp2)); //subject name: Query name: subject start: subject end : Query start :  Query end: Query length : subject length
                    else { // same subject found for another contig
                        int a = std::distance(entry_blast.begin(), entries_it);
                        blast_temp = entry_blast[temp[2]];
                        string contig1 = blast_temp[0];
                        sub_st = atoi(blast_temp[1].c_str());
                        sub_en = atoi(blast_temp[2].c_str());
                        sub_st1 = atoi(temp[13].c_str());
                        sub_en1 = atoi(temp[14].c_str());

                        if (contig1.substr(0, contig1.size() - 3) == temp[0].substr(0, temp[0].size() - 3)) { //BOTH CONTIGS FROM SAME CLUSTER AN D HAVE SAME NAME, ignore them
                            entry_blast.insert(make_pair(temp[2], temp2)); //subject name: Query name: subject start: subject end 
                            check = true;
                            continue;
                        } else if ((sub_st1 < sub_st && sub_en1 < sub_st) || (sub_st1 > sub_en && sub_en1 > sub_en)) { // map on different places on subject
                          //  if()
                            blast_analysis << temp[2] << "\t" << temp[3] << "\t" << contig1 << "\t" << blast_temp[5] << "\t" << sub_st << "\t" << sub_en << "\t" << temp[0] << "\t" << temp[1] << "\t" << temp[13] << "\t" << temp[14] << endl;
                            entry_blast.erase(entries_it);
                        } else {
                            entry_blast.insert(make_pair(temp[2], temp2));
                        }
                    }
                    id++;
                    check = false;
                    continue;
                }
                temp2.clear();
                blast_temp.clear();
            } else {
                check = true;
                id = 0;
            }
        }
        exit(0);
        /*ALl possible partial contigs*/
     /*   
        map<string, vector<string> >::iterator it;
        for (it = entry_blast.begin(); it != entry_blast.end(); it++) {
            blast_temp = it->second;
            if (temp[1] < temp[2]) {
                sub_st = atoi(blast_temp[1].c_str());
                sub_en = atoi(blast_temp[2].c_str());
            } else {
                sub_st = atoi(blast_temp[2].c_str());
                sub_en = atoi(blast_temp[1].c_str());
            } //reverse strand
            q_st = atoi(blast_temp[3].c_str());
            q_end = atoi(blast_temp[4].c_str());
            q_len = atoi(blast_temp[5].c_str());
            sub_len = atoi(blast_temp[6].c_str());
            query_cov = float((q_len / sub_len))*100;
            if (query_cov <= 90) { //if((q_st - sub_st >=10) ||(sub_en - q_end >=10)){ // to get partial 
                blast_analysis << blast_temp[0] << "\t" << blast_temp[3] << "\t" << blast_temp[4] << "\t" << blast_temp[1] << "\t" << blast_temp[2] << endl;
            } else
                continue;
        } */
        //entry_blast.clear();
        string contig;
        /*ALl possible false chimeras contigs*/
       // blastn.close();
        //blastn.open(blastn_file.c_str());
        while (!(getline(blastn, entry).eof())) {
            if (entry[0] == '#' && check) {
                entry_blast.clear();
                id = 0;
                //check = false;
                continue;
            } 
            else if (entry[0] == '#' && id < 2)  { // for 1 entry
                id = 0;
                check = true;
                continue;
            } 
            else if (id < 2) {
                vector<string> temp3 = split(entry, "\t");
                query_cov = atoi(temp3[4].c_str());
                q_len = atoi(temp3[1].c_str());
                contig = temp3[0];

                vector<string> temp2;
                temp2.push_back(temp3[2]);
                temp2.push_back(temp3[11]);
                temp2.push_back(temp3[12]); // temp2.push_back(temp[11]); temp2.push_back(temp[12]);
                
                if (query_cov < query_cov_TH) {
                    entry_blast.insert(make_pair(temp3[2], temp2)); //subject name: Query name: subject start: subject end     id++;                       
                }
                check = false;
                id++; // continue;
            } else if (entry_blast.size() == 2 && !check) { // two blast hits to be false chimeras
                it = entry_blast.begin();
                blast_temp = it->second;
                it++;
                blast_temp2 = it->second;
                //               blast_temp = entry_blast[0];
                //             blast_temp2 = entry_blast[1];
                q_st = atoi(blast_temp[1].c_str());
                q_end = atoi(blast_temp[2].c_str());
                q_st1 = atoi(blast_temp2[1].c_str());
                q_end1 = atoi(blast_temp2[2].c_str());

                if ((blast_temp[0] != blast_temp2[0]) && ((q_st1 < q_st && q_end1 < q_end) || (q_st1 > q_st && q_end1 > q_end))) { // map on different places on subject
                    blast_analysis << contig << "\t" << q_len << "\t" << blast_temp[0] << "\t" << q_st << "\t" << q_end << "\t" << blast_temp2[0] << "\t" << q_st1 << "\t" << q_end1 << endl;
                }

                check = true;
                id = 0;
                blast_temp.clear();
                blast_temp2.clear();
                entry_blast.clear();
            }
            else {
                check = true;
            }
            //                            
        }
        blast_analysis.close();
        blastn.close();
    }

    void repeat_contig(string fasta_file, string path_file, string blast_file) {

        //  map<string, map<string, string> > blast_results;

        typedef map<string, string> Map1; // inner map
        typedef map<string, Map1> Map2; //  middle map
        typedef map<string, Map2> Map3; //  outer map
        Map3 blast_results;

        string entry1, entry2, entry3, id_first = "", id, path;
        vector<string> path1, path2;
        vector<string> temp, temp2, path_map, check_path;
        map<string, string> fasta_header, id_path;
        float count = 0, count1 = 0, match = 0;
        bool found = false;
        typedef multimap<string, string > mmap;
        mmap merged;
        typedef pair<string, string> id1_id2;
        ifstream fastafile, pathfile, blastfile;
        blastfile.open(blast_file.c_str());
        while (!(getline(blastfile, entry3)).eof()) {
            if (entry3 [0] == '#') {
                count1 = 0;
                continue;
            } else if (count1 <= 3) {
                temp2 = split(entry3, "\t");
                // blast_results.insert(make_pair(temp2[0], Map1())); //nested map with multiple valued per key
                blast_results[temp2[0]][temp2[1]].insert(make_pair(temp2[6], temp2[9])); //ID1 ID2 ID1_start ID2_start
                // blast_results[temp2[0]].insert(make_pair(temp2[1], temp2[9]));
                count1++;
            }
        }

        cout << blast_results.size() << endl;
        blastfile.close();
        fastafile.open(fasta_file.c_str());
        pathfile.open(path_file.c_str());
        while (!(getline(pathfile, entry1)).eof()) { //read path file

            if (entry1[0] == '>') {
                std::size_t found = entry1.find("path");
                temp = split(entry1, " ");
                string id_fasta = temp[0].substr(1, temp[0].length() - 4);
                string path_fasta = entry1.substr(found, entry1.length() - 1);
                fasta_header.insert(make_pair(id_fasta, path_fasta));
                //    [id_fasta] = path_fasta;

            }
        }
        pathfile.close();

        while (!(getline(fastafile, entry2).eof())) {
            if (entry2[0] == '>') {
                id = entry2.substr(1, entry2.size() - 1);
                if (count == 0) {
                    id_first = id;
                    count = 1;
                } else {
                    if (id.substr(0, id.size() - 3) == id_first.substr(0, id_first.size() - 3)) {
                        check_path.push_back(id_first);
                        id_first = id;
                    } else { //collect genes of same cluster 
                        // same_contig.push_back(first_entry);
                        check_path.push_back(id_first);

                        for (int a = 0; a < check_path.size(); a++) {
                            map<string, string> ::iterator itt = fasta_header.find(check_path[a]);
                            //find(check_path [a]);
                            if (itt != fasta_header.end()) { //if found add whole entry in new vector to process further
                                vector <string> m = split(itt->second, "]");
                                path = m[1].substr(5, m[1].size() - 9); // extract last [] from path nodes
                                id_path[check_path[a]] = path;
                                found = true;
                            } else if (found == true) break;
                            else continue;
                        }
                        map<string, string>::iterator itr_path;
                        map<string, string> merge_id;
                        string id1, id2, seq1, seq2, merge_seq;
                        float match_th = 0.0, j = 0, i = 0, start_id1 = 0, start_id2 = 0;

                        //check for path
                        for (itr_path = id_path.begin(); itr_path != id_path.end(); itr_path++) {
                            if (i == id_path.size() - 1) break;
                            cout << itr_path ->first << endl;
                            id1 = itr_path->first;
                            path1 = split(itr_path->second, ",");
                            map<string, string>::iterator itr_path2 = id_path.begin();
                            i++;
                            j = i;

                            for (std::advance(itr_path2, j); itr_path2 != id_path.end(); itr_path2++) {
                                j++;
                                path2 = split(itr_path2->second, ",");
                                id2 = itr_path2->first;

                                match_th = (path2.size() / 2.0);
                                cout << path2.size() << "  " << match_th << endl;

                                //acc to observation path node of min path ID is considered as standard
                                if (path2.size() > path1.size()) {
                                    match_th = ceil(path1.size() / 2.0);
                                    for (int c = 0; c < path1.size(); c++) {
                                        ptrdiff_t it = std::find(path2.begin(), path2.end(), path1[c]) - path2.begin();
                                        if (it != path2.size()) {
                                            match++;
                                        }
                                    }
                                } else {
                                    match_th = ceil(path2.size() / 2.0);
                                    for (int c = 0; c < path2.size(); c++) {
                                        ptrdiff_t it = std::find(path1.begin(), path1.end(), path2[c]) - path1.begin();
                                        if (it != path1.size()) {
                                            match++;
                                        }
                                    }
                                }
                                multimap<string, string > ::iterator itt_id;
                                if (match >= match_th) {
                                    //   multimap<string, string > ::iterator itt_id1 = merged.find(id1);
                                    if (!(merged.empty())) {
                                        for (itt_id = merged.begin(); itt_id != merged.end(); ++itt_id)
                                            if (itt_id->first == id1) {
                                                //  if (itt_id1!= merged.end()) { //merge id1 and id2 // dint merged before
                                                //merged[id1] = id2;
                                                merged.insert(id1_id2(itt_id->first, id2));
                                                break;
                                                //                                        merged.push_back(id2);
                                                //                                        seq1 = extract_fasta(id1, fasta_file);
                                                //                                        seq2 = extract_fasta(id2, fasta_file);
                                                //                                        merge_seq = align_merge(merge_id[id1], seq2, start_id1, start_id2);
                                                //                                        merge_id[id1] = merge_seq;
                                                //                                        merge_id[id2] = merge_seq;

                                            } else if (itt_id->second == id1) { //merge id1 and id2 // merged before so use merged sequence for further merging
                                                //                                        merged.push_back(id1);
                                                merged.insert(id1_id2(itt_id->first, id2));
                                                break;
                                                //                                        seq1 = extract_fasta(id1, fasta_file);
                                                //                                        seq2 = extract_fasta(id2, fasta_file);
                                                //                                        merge_seq = align_merge(seq1, merge_id[id2], start_id1, start_id2);
                                                //                                        merge_id[id1] = merge_seq;
                                                //                                        merge_id[id2] = merge_seq;

                                            } else if (itt_id->second == id2) { //merge id1 and id2 // dint merged before
                                                merged.insert(id1_id2(itt_id->first, id1));
                                                break;
                                                //                                        merged.push_back(id1);
                                                //                                        merged.push_back(id2);
                                                //                                        seq1 = extract_fasta(id1, fasta_file);
                                                //                                        seq2 = extract_fasta(id2, fasta_file);
                                                //                                        merge_seq = align_merge(seq1, seq2, start_id1, start_id2);
                                                //                                        merge_id[id1] = merge_seq;
                                                //                                        merge_id[id2] = merge_seq;
                                            } else {
                                                merged.insert(id1_id2(id1, id2));
                                                break;

                                            }

                                    } else {
                                        merged.insert(id1_id2(id1, id2));

                                    }
                                    match = 0;
                                    //        cout << id1 << "  " << id2 << endl << merge_seq << endl;
                                }
                                //else no path node match
                            }
                        }
                        //search blastn data for start position of alignment   
                        map<string, Map2>::iterator it_out = blast_results.find(id1);
                        map<string, Map1>::iterator it_mid;
                        map<string, string>::iterator it_in;


                        for (it_mid = (*it_out).second.begin(); it_mid != (*it_out).second.end(); it_mid++) {
                            //output the 1st and 2nd elements of Map1
                            if ((*it_mid).first == id2) {
                                it_in = (*it_mid).second.begin();
                                string pos = (*it_in).first;
                                string pos1 = (*it_in).second;
                                start_id1 = atoi(pos.c_str());
                                start_id2 = atoi(pos1.c_str());
                                break;
                            } else
                                start_id1 = 0;
                            start_id2 = 0;
                        }
                        for (mmap::iterator it = merged.begin(); it != merged.end(); ++it) {

                            std::cout << it->first << " " << it->second << endl;
                        }
                        cout << "hello" << endl;
                        merged.clear();
                        merge_id.clear();
                        check_path.clear();
                        found = false;
                        id_first = id;
                    }

                }
            }
        }

    }

    string align_merge(string str1, string str2, int start_id1, int start_id2) {
        int discard_seq_TH = 20;
        int n, m, pos_MI = 1, pos_RI;
        int k = 1, consecutive_mm = 0;
        int counter = 0;
        int max_mismatch = 20, mismatch = 0;
        string overlap = "", maxOverlap = "", final, MI, RI, truncated_RI, truncated_MI, before_aln = "";
        int overlapSize = 0;
        int i = start_id1, j = start_id2;
        if ((start_id1 > 1) && (start_id1 > start_id2)) before_aln = str1.substr(0, start_id1 - 1);
        else if ((start_id2 > 1) && (start_id2 > start_id1)) before_aln = str2.substr(0, start_id2 - 1);


        char base_str1, base_str2;
        while (i < str1.size()) {
            while (j < str2.size()) {
                base_str1 = str1[i];
                base_str2 = str2[j];
                if (base_str1 == base_str2) {
                    overlap = overlap + base_str1;
                    pos_MI = j;
                    i++;
                    j++;
                    consecutive_mm = 0;
                    goto out_of_j;
                } else if ((j < str2.length() - 4) && (i < str1.length() - 4) && (str1[i] == str2[j + 1]) && (str1[i + 1] == str2[j + 2]) && (str1[i + 2] == str2[j + 3])) { // indel check, 3 bases check for conformation
                    overlap = overlap + base_str2; //insertion in string2 
                    j++;
                    pos_MI = j;
                    goto out_of_j;
                } else if ((j < str2.length() - 4) && (i < str1.length() - 4) && (str1[i + 1] == str2[j]) && (str1[i + 2] == str2[j + 1]) && (str1[i + 3] == str2[j + 2])) { // indel check
                    overlap = overlap + base_str1; //insertion in string1 
                    i++;
                    pos_MI = j;
                    goto out_of_j;
                } else {

                    char code = base_code(base_str1, base_str2);
                    mismatch++;
                    if (mismatch <= max_mismatch && consecutive_mm <= 5) {
                        overlap = overlap + code;
                        consecutive_mm++;
                        i++;
                        j++;
                        pos_MI = j;

                        goto out_of_j;
                    } else {
                        if (str1.size() > str2.size()) truncated_MI = str1.substr(pos_MI + 1, str1.size() - 1);
                        else truncated_MI = str2.substr(pos_MI + 1, str2.size() - 1);
                        pos_MI = j; // size of overlap
                        //                    pos_RI = i;
                        //                    RI = str1.substr(0, str1_st);
                        //                    MI = str2.substr(maxOverlap.size(), str2.size() - 1);
                        //                    truncated_RI = str1.substr(pos_RI + 1, str1.size() - 1); //end of contig1 left from overlap
                        //                    truncated_MI = str2.substr(pos_MI + 1, str2.size() - 1);
                        //                    std::transform(maxOverlap.begin(), maxOverlap.end(), maxOverlap.begin(), ::tolower); //end of contig2 left from overlap
                        final = before_aln + overlap + truncated_MI; //RI + maxOverlap + MI;
                        goto out_of_i;
                    }
                }
            }
out_of_j:
            ;
        }
out_of_i:
        ;
        if (str1.size() > str2.size()) truncated_MI = str1.substr(pos_MI + 1, str1.size() - 1);
        else truncated_MI = str2.substr(pos_MI + 1, str2.size() - 1);
        final = before_aln + overlap + truncated_MI;
        overlap.clear();
        if (final != "") {
            //    cout << str1 << endl;
            //    cout << str2 << endl;
            //   cout << final << endl;
        } else {
            cout << "no overlap found" << endl;
            return final;
        }
    }

    char base_code(char base1, char base2) {
        char code;
        if ((base1 == 'G' && base2 == 'A') || (base1 == 'A' && base2 == 'G')) code = 'R'; // puRine
        else if ((base1 == 'T' && base2 == 'C') || (base1 == 'C' && base2 == 'T')) code = 'Y'; //pYrimidine
        else if ((base1 == 'A' && base2 == 'C') || (base1 == 'C' && base2 == 'A')) code = 'M'; //aMino
        else if ((base1 == 'G' && base2 == 'T') || (base1 == 'T' && base2 == 'G')) code = 'K'; //Keto
        else if ((base1 == 'G' && base2 == 'C') || (base1 == 'C' && base2 == 'G')) code = 'S'; //Strong interaction (3 H bonds)
        else if ((base1 == 'A' && base2 == 'T') || (base1 == 'T' && base2 == 'A')) code = 'W'; //Weak interaction (2 H bonds)
        else code = 'X';
        return code;


    }

    void remove_file(const string& file) {
        remove(file.c_str());
    }

    map<string, map<int, int> > coverage_data(string cov_input) // written for coverage file obtained from samtools depth
    {
        string entry, cur_ID, prev_ID;
        int count = 0;
       map<string, map<int, int> > coverage;
        map<int, int> cov;
        vector<string> temp1;
        ifstream cov_file(cov_input.c_str());
        while (!(getline(cov_file, entry).eof())) {
            temp1 = split(entry, "\t"); // ID, source, feature, start, stop, score, strand/feature, attributes
            //str_split(entry, temp1, delimiter);
            coverage[temp1[0]][atoi(temp1[1].c_str())] = atoi(temp1[2].c_str());   
        }
        cov_file.close();
        return coverage;
    }

    map<string, vector<int> > coverage_data_updated(string cov_input) // written for coverage file obtained from samtools depth
    {
        string entry, cur_ID, prev_ID;
        int count = 0;
        map <string, vector<int> > coverage;
        vector<int> cov;
        vector<string> temp1;
        ifstream cov_file(cov_input.c_str());
        while (!(getline(cov_file, entry).eof())) {
            temp1 = split(entry, "\t"); // ID, source, feature, start, stop, score, strand/feature, attributes

            if (count == 0) {
                prev_ID = temp1[0];
                cov.push_back(atoi(temp1[2].c_str())); //take coverage only, as position should be same as index number
                count++;
            } else {
                cur_ID = temp1[0];
                if (cur_ID == prev_ID) {
                    cov.push_back(atoi(temp1[2].c_str())); //take coverage only, as position should be same as index number
                    prev_ID = cur_ID;

                } else {
                    coverage[prev_ID] = cov;
                    prev_ID = cur_ID;
                    cov.clear();
                    cov.push_back(atoi(temp1[2].c_str())); //add cov of new contig
                }
            }
        }
        cov_file.close();
        return coverage;
    }

    int average_insertsize(string bam_file) { // average insert size of first contig 154 for arabidopsis 1655112
        BamReader reader;
        bool pad;
        float count = 0.0, insert_size = 0.0, avg_is = 0.0;
        int prev_contigID = 0;
        // Opens a bam file
        if (!reader.Open(bam_file)) {
            cerr << "Could not open BAM file." << endl;
            exit(0);
        }
        // RefData read_id, mate_id;
        vector<RefData> ref;
        ref = reader.GetReferenceData();
        BamTools::BamAlignment al;
        while (reader.GetNextAlignment(al)) {
            read_length = al.Length;
            if (al.RefID == prev_contigID) { // to avoid unaligned reads which gives -1 RefID
                if (al.RefID == al.MateRefID) {
                    insert_size = insert_size + abs(al.InsertSize);
                    count++;
                }
                prev_contigID = al.RefID;
            } else {
                avg_is = float(insert_size / count);
                prev_contigID = al.RefID;
                break;
            }
        }

        reader.Close();
        return avg_is;
    }

    map<string, list<int> > cufflink_data(string cufflink_input) {
        string entry, cur_ID, prev_ID;
        int count = 0;
        map <string, list<int> > split_site;
        list<int> exon;
        vector<string> temp1;
        ifstream cufflink_file(cufflink_input.c_str());
        while (!(getline(cufflink_file, entry).eof())) {
            temp1 = split(entry, "\t"); // ID, source, feature, start, stop, score, strand/feature, attributes

            if (count == 0) {
                prev_ID = temp1[0];
                if (temp1[2] == "exon") // ignore transcript positions
                {
                    exon.push_back(atoi(temp1[3].c_str()));
                    exon.push_back(atoi(temp1[4].c_str()));
                }
                count++;
            } else {
                cur_ID = temp1[0];
                if (cur_ID == prev_ID) {
                    if (temp1[2] == "exon") // ignore transcript positions
                    {
                        exon.push_back(atoi(temp1[3].c_str()));
                        exon.push_back(atoi(temp1[4].c_str()));
                    }
                    prev_ID = cur_ID;

                } else {
                    split_site[prev_ID] = exon;
                    exon.clear();
                    prev_ID = cur_ID;
                    if (temp1[2] == "exon") { // add data for new contig
                        exon.push_back(atoi(temp1[3].c_str()));
                        exon.push_back(atoi(temp1[4].c_str()));
                    }

                }
            }
        }
        cufflink_file.close();
        return split_site;
    }

    struct length {

        bool operator()(const string& a, const string& b) {
            return a.size() < b.size();
        }
    };

string consensus_seq_CAP3(std::vector<string>softclip_list, bool toleft, string PATH) {

    string SCs_file = PATH + "/SCs_file.fa";
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
         //   cout << endl << "more than 1 consensus sequence found for extension using softclips " << endl;
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
        } else { // if CAP3 doesn't generate any consensus sequennce (for shorter sequences), call consensus_seq function to generate consensus seq manually
            cons_seq = consensus_seq(softclip_list, toleft);
        }


remove(CAP3_file.c_str());
remove(SCs_file.c_str());
return cons_seq;

}
      string consensus_seq(std::vector<string>softclip_list, bool left) {

        float cov_a = 0, consec_a = 0; // for poly A tails
        float cov_c = 0;
        float cov_g = 0;
        float cov_t = 0;
        int size = 0;
        float count = 0;
        string final_cons_sequence;
        int count_amb_code = 0;
        int max_amb_code_allowed = 3;
        int full_A_tail = 90;
        int partial_A_tail = 40;
        int max_cov = 60;
        string str = "";
        char code;
        sort(softclip_list.begin(), softclip_list.end(), length());
        string s = softclip_list.back();
        //string s = softclip_list[softclip_list.size()-1];
        size = s.size() - 1;
        vector <char> con_seq; // [size-1] = "";

        // cout << size << flag << endl;
        //  if (size > 1) {
        if (left) { //start

            int t = size;
            int clip_size = 1;
            while (t > 0) {
                for (int z = 0; z < softclip_list.size(); z++) {
                    str = softclip_list[z];
                    int itt = str.size(); // start from end of sc or start of read
                    if (itt < clip_size)
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
                // cout<<cov_a<<" "<<cov_c<<" "<<cov_g<<" "<<cov_t<<endl;
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
                     /*   if (cov_c == cov_a && cov_c > cov_g && cov_c > cov_t) {
                            code = base_code('C', 'A');
                            count_amb_code++;
                            con_seq.push_back(code);
                        } else if (cov_g == cov_a && cov_g > cov_c && cov_g > cov_t) {
                            code = base_code('G', 'A');
                            count_amb_code++;
                            con_seq.push_back(code);
                        } else if (cov_t == cov_a && cov_t > cov_c && cov_t > cov_g) {
                            code = base_code('T', 'A');
                            count_amb_code++;
                            con_seq.push_back(code);
                        } else if (cov_g == cov_c && cov_g > cov_a && cov_c > cov_t) {
                            code = base_code('G', 'C');
                            count_amb_code++;
                            con_seq.push_back(code);
                        }//
                        else if (cov_g == cov_t && cov_g > cov_c && cov_g > cov_a) {
                            code = base_code('G', 'T');
                            count_amb_code++;
                            con_seq.push_back(code);
                        } else if (cov_t == cov_c && cov_t > cov_a && cov_t > cov_g) {
                            code = base_code('C', 'T');
                            count_amb_code++;
                            con_seq.push_back(code);
                        } else { */
                            code = 'N';
                            count_amb_code++;
                            con_seq.push_back(code);

                      //  }
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
                        break; //don't allow 3 base mismatches or N
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
                // cout<<cov_a<<" "<<cov_c<<" "<<cov_g<<" "<<cov_t<<endl;
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
                        /*if (cov_c == cov_a && cov_c > cov_g && cov_c > cov_t) {
                            code = base_code('C', 'A');
                            count_amb_code++;
                            con_seq.push_back(code);
                        }//
                        else if (cov_g == cov_a && cov_g > cov_c && cov_g > cov_t) {
                            code = base_code('G', 'A');
                            count_amb_code++;
                            con_seq.push_back(code);
                        } else if (cov_t == cov_a && cov_t > cov_c && cov_t > cov_g) {
                            code = base_code('T', 'A');
                            count_amb_code++;
                            con_seq.push_back(code);
                        } else if (cov_g == cov_c && cov_g > cov_a && cov_c > cov_t) {
                            code = base_code('G', 'C');
                            count_amb_code++;
                            con_seq.push_back(code);
                        }//
                        else if (cov_g == cov_t && cov_g > cov_c && cov_g > cov_a) {
                            code = base_code('G', 'T');
                            count_amb_code++;
                            con_seq.push_back(code);
                        } else if (cov_t == cov_c && cov_t > cov_a && cov_t > cov_g) {
                            count_amb_code++;
                            code = base_code('C', 'T');
                            con_seq.push_back(code);
                        } else { */
                            count_amb_code++;
                            code = 'N';
                            con_seq.push_back(code);

                       // }
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
        if (final_cons_sequence.length() < 5)
            final_cons_sequence = "";
        //  }
        return final_cons_sequence;
    }
    void simulated_errors(string fasta_file, string sim_err_fastafile, string sim_err_log) {
        utils utils;
        string entry, ID, fasta_seq, first_seq, sec_seq, false_chimera, id1, id2, ID1;
        list<string> all_IDs, simulated_IDs;
        ofstream err_fasta;
        err_fasta.open(sim_err_fastafile.c_str());
        ofstream err_log;
        err_log.open(sim_err_log.c_str());
        ifstream fastafile;
        fastafile.open(fasta_file.c_str());
        int j = rand(), chimera_pair = 1, chimera_count = 1, broken_count = 1;
        while (!(getline(fastafile, entry).eof())) {
            if (entry[0] == '>') {
                all_IDs.push_back(entry.substr(1, entry.size() - 1));
            }
        }
        fastafile.close();
        std::list<std::string>::iterator it = all_IDs.begin();
        int i = 0;
        while (i <= 999) {
            float false_chimeric_pos1 = rand() % 20 + 40; // 40-60% of the length of contig1 to be merged with contig2
            float false_chimeric_pos2 = rand() % 20 + 50; // 50-70% of the length of contig1 to be merged with contig2
            int random_ID = rand() % all_IDs.size() + 1;
            std::advance(it, random_ID);
            ID = *it;
            fasta_seq = utils.extract_fasta(ID, fasta_file);
            int size = fasta_seq.size();
            if (chimera_count <= 500) {
                if (chimera_pair == 1) {
                    false_chimeric_pos1 = (false_chimeric_pos1 * fasta_seq.size()) / 100;
                    first_seq = fasta_seq.substr(0, false_chimeric_pos1);
                    id1 = ID;
                    chimera_pair++;
                } else if (chimera_pair == 2) {
                    false_chimeric_pos2 = (false_chimeric_pos2 * fasta_seq.size()) / 100;
                    sec_seq = fasta_seq.substr(false_chimeric_pos2, fasta_seq.size()); // take the later half to merge because read coverage in the first half is usually low so it will not give right pattern in the middle
                    id2 = ID;
                    chimera_pair++;

                } else {
                    chimera_pair = 1;
                    false_chimera = first_seq + sec_seq;
                    err_fasta << ">" << id1 << "_" << id2 << endl << false_chimera << endl;
                    err_log << id1 << "\t" << id2 << "\t" << "chimera" << endl;
                    chimera_count++;
                    i++;
                    simulated_IDs.push_back(id1);
                    simulated_IDs.push_back(id2);
                }

            } else if (broken_count <= 500 && size >= 1000) {
                float false_broken = rand() % 20 + 40; // 40-60% of the length of contig1 to be broken
                false_broken = (false_broken * fasta_seq.size()) / 100;
                err_fasta << ">" << ID << "_X" << endl << fasta_seq.substr(0, false_broken) << endl;
                err_fasta << ">" << ID << "_Y" << endl << fasta_seq.substr(false_broken, fasta_seq.size()) << endl;
                err_log << ID << "\t" << false_broken << "\t" << "fragmented" << endl;
                broken_count++;
                i++;
                simulated_IDs.push_back(ID);
            }

        }

        std::list<std::string>::iterator it_find;
        for (it = all_IDs.begin(); it != all_IDs.end(); it++) {
            it_find = find(simulated_IDs.begin(), simulated_IDs.end(), *it);
            if (it_find == simulated_IDs.end()) {
                fasta_seq = utils.extract_fasta(*it, fasta_file);
                err_fasta << ">" << *it << endl << fasta_seq << endl;
            }
        }
        err_fasta.close();
        err_log.close();
    }

    void fastq_dataset(string fasta_in, string bam_in, string fastq_first_out, string fastq_sec_out) {
        BamReader reader;
        size_t found = fasta_in.find_last_of("/\\");
        string set_path = fasta_in.substr(0, found);
        string bam_out = set_path + "/new.bam";
        vector<int> clipsize, read_pos, gen_pos;
        bool pad;
        int FLAG_READ_UNMAPPED = 4;
        int FLAG_MATE_UNMAPPED = 8;
        int count = 0;
        set <string> reads_to_keep;

        string ID_temp, entry_fasta;
        vector<string> ID;
        ifstream fastafile;
        fastafile.open(fasta_in.c_str());
        while (!(getline(fastafile, entry_fasta).eof())) {
            if (entry_fasta[0] == '>') {
                temp = split(entry_fasta, " ");
                ID_temp = entry_fasta;
                ID.push_back(temp[0].substr(1, temp[0].size() - 1));
            }
        }



        const SamHeader header = reader.GetHeader(); // returns header data object
        // Opens a bam file
        if (!reader.Open(bam_in)) {
            cerr << "Could not open BAM file." << endl;
            exit(0);
        }
        // RefData read_id, mate_id;
        vector<RefData> ref;
        BamTools::RefData ref_contig, mate_contig;
        ref = reader.GetReferenceData();

        BamWriter writer;
        if (!writer.Open(bam_out, header, ref)) {
            cerr << "Could not open output BAM file" << endl;
            exit(0);
            // return;
        }

        BamTools::BamAlignment al;
        while (reader.GetNextAlignment(al)) {
            if (al.RefID != -1) {
                //      cigar = al.CigarData;
                string read_name = al.Name;
                ref_contig = ref.at(al.RefID);
                int len_read_contig = ref_contig.RefLength;

                BamTools::RefData mate_contig = ref.at(al.MateRefID);
                int len_mate_contig = mate_contig.RefLength;
                if (reads_to_keep.find(read_name) != reads_to_keep.end()) // read already added, add its mate too(mapped or unmapped both)
                {
                    writer.SaveAlignment(al);
                    count++;
                    reads_to_keep.insert(read_name);
                } else { //not already added? 
                    if (std::find(ID.begin(), ID.end(), ref_contig.RefName) != ID.end()) { //current contig is not being merged

                        writer.SaveAlignment(al);
                        reads_to_keep.insert(read_name);
                        count++;
                    } else
                        continue;
                }
            }
        }
        reader.Close();
        writer.Close();
        //emove(bam_in.c_str());
        filterSamToFastq ff;
        ff.filterFastq(bam_out, fastq_first_out, fastq_sec_out);
    }

    void str_split(const string& str, vector<string>& tokens, const string& delimiters) {
        // Skip delimiters at beginning.
        string::size_type last_pos = str.find_first_not_of(delimiters, 0);
        // Find first "non-delimiter".
        string::size_type pos = str.find_first_of(delimiters, last_pos);

        while (string::npos != pos || string::npos != last_pos) {
            // Found a token, add it to the vector.
            tokens.push_back(str.substr(last_pos, pos - last_pos));
            // Skip delimiters.  Note the "not_of"
            last_pos = str.find_first_not_of(delimiters, pos);
            // Find next "non-delimiter"
            pos = str.find_first_of(delimiters, last_pos);
        }
    }

    void findPSC_byBLAST(string blastn_file, string potential_PSC_file) { // not whole assembly blast instead blast only identified contigs 

        ifstream blastn;

        string entry_blast, subject_ID, ID, query_ID, subject_ID_prev;

        float aln_len, query_len, p_identity, accuracy, e_value, subject_len;

        bool check = true, global_transcript = false;

        vector<string> temp;

        list<string> subject_list;

        float count = 0.0, total_merged = 0.0;

        blastn.open(blastn_file.c_str());

        ofstream out_file;

        out_file.open(potential_PSC_file.c_str());

        while (!(getline(blastn, entry_blast).eof())) {

            if (entry_blast[0] == '#') {

                if (subject_list.size() > 1) {

                    out_file << ID << "\t";

                    cout << ID << "\t";

                    for (list<string>::iterator it = subject_list.begin(); it != subject_list.end(); it++)

                        out_file << *it << "\t";

                    out_file << endl;
                }
                subject_list.clear();

                check = true;

                global_transcript = false;

                continue;

            } else if (check == true) {

                temp.clear();

                str_split(entry_blast, temp, delimiter);

                ID = temp[0];

                subject_ID = temp[2];

                subject_len = atof(temp[3].c_str());
                e_value = atof(temp[15].c_str());

                //                p_identity = atof(temp[6].c_str()); 
                //
                //                matched_bases = atoi(temp[6].c_str()); 

                aln_len = atof(temp[7].c_str());

                query_len = atof(temp[1].c_str());

                float coverage = atof(temp[4].c_str()); //(aln_len / query_len)*100; 

                accuracy = atof(temp[6].c_str()); //(matched_bases/aln_len) *100 ; 

                if (coverage >= 90 && accuracy >= 75 && global_transcript == false && fabs(subject_len - query_len) <= 10) {

                    global_transcript = true;

                    check = true;

                    count++;

                    subject_list.push_back(subject_ID);
                    subject_ID_prev = subject_ID;

                } else if (e_value <= 0.001 && accuracy >= 70 && global_transcript == true) {
                    if (subject_ID != subject_ID_prev) {

                        subject_ID_prev = subject_ID;

                        subject_list.push_back(subject_ID);

                        check = true;

                        count++;
                    }
                } else {
                    global_transcript = false;

                    check = false;
                }

            }

        } //end reading blastn file 

        blastn.close();

        out_file.close();

    }

    void heterozygosity_count(string vcf_file, string heterozygosity_count_file) { // not whole assembly blast instead blast only identified contigs 

        ifstream vcf;
        vcf.open(vcf_file.c_str());
        ofstream count_file;
        count_file.open(heterozygosity_count_file.c_str());
        string ID, ID_prev, entry_vcf;
        int count = 0, heterozygous_SNP_count = 0;
        vector<string> temp;
        
        while (!(getline(vcf, entry_vcf).eof())) {
            temp.clear();
            str_split(entry_vcf, temp, delimiter);

            if (count == 0) {
                ID_prev = temp[0];
            } else {
                ID = temp[0];
                if (ID == ID_prev) {
                    heterozygous_SNP_count++;
                    ID_prev = ID;
                } else {
                    heterozygous_SNP_count++;
                    count_file << ID_prev << "\t" << heterozygous_SNP_count << endl;
                    ID_prev = ID;
                    heterozygous_SNP_count = 0;
                }
            }
            count++;

        }
        vcf.close();
        count_file.close();
    }

    int findNthOccur(string str, char ch, int N) //start from 0th position
    {
        int occur = 0;
        vector<int> occur_ind;

        // Loop to find the Nth 
        // occurence of the character 
        for (int i = 0; i < str.length(); i++) {
            if (str[i] == ch) {
                // occur += 1; 
                occur_ind.push_back(i);
            }
        }
        if (!(occur_ind.empty())) {
            if (N < 0)
                occur = occur_ind.end()[N];
            else
                occur = occur_ind.begin()[N];
            return occur;
        }
        return -1;
    }

    vector<string> str_split2(string s, string delimiter) {
        vector<string> out;
        //std::string s = "scott>=tiger>=mushroom";
        //std::string delimiter = ">=";

        size_t pos = 0;
        std::string token;
        while ((pos = s.find(delimiter)) != std::string::npos) {
            token = s.substr(0, pos);
            out.push_back(token); // << std::endl;
            s.erase(0, pos + delimiter.length());
        }
        //std::cout << s << std::endl;
        out.push_back(s);
        return out;
    }

 void format_finalAssembly(string in_assembly, string out_assembly , string initial_IDs_map) {

        //avg_IS  = 120;
        avg_IS  = 5;
        std::size_t found = out_assembly.find_last_of("/\\");
        size_t find_ROAST;
        string summary_file = out_assembly.substr(0, found) + "/summary.txt";

                        
        ifstream in_file;
        in_file.open(in_assembly.c_str());
        
        ifstream init_map;
        init_map.open(initial_IDs_map.c_str());

        ofstream out_file;
        out_file.open(out_assembly.c_str());
        
        ofstream summary;
        summary.open(summary_file.c_str());

        ofstream del_seq;

        if (ignore_short_seq > 0) {
            string remove_seq = out_assembly.substr(0, found) + "/removed_seq.fasta";
            del_seq.open(remove_seq.c_str());
        }

        string entry, ID, ID2, fasta_seq, final_fasta_seq;
        map<string, list<string> >linked_contigs;
        list<string> temp;
        list <string> temp_ID;
        vector<string> IDs, temp_ID2;
        map<string, string> init_IDs_map;
        stringstream ss;
        int ind = 0, ind2 = 0, contig_num = 1;

        while (!getline(in_file, entry).eof()) { //save all IDs of fastafile
            if (entry[0] == '>') {
                ID = entry.substr(1, entry.size() - 1);
                //  cout << ID <<endl;
                IDs.push_back(ID);
            }
        }
        in_file.close();

    while (!getline(init_map, entry).eof()) { //save data from initial IDs map
            
            str_split(entry, temp_ID2, "\t");
            init_IDs_map[temp_ID2[0]] = temp_ID2[1];
            temp_ID2.clear();
        }
        init_map.close();
        temp_ID2.clear();
        // cout << "got all IDs" << endl;
        vector<string> IDs2(IDs);

        for (std::vector<std::string>::iterator it = IDs.begin(); it != IDs.end(); ++it) {
            ID = *it;
            // cout << ID <<endl;
            if (ID.find(ROASTcap3_left_tag) == string::npos && ID.find(ROASTcap3_right_tag) == string::npos) { // if cap3 not found , get unique ID
                
                for (std::vector<std::string>::iterator it2 = IDs2.begin(); it2 != IDs2.end(); it2++) {
                    
                    ID2 = *it2;
                    int find_cp3;
                    find_cp3 = ID2.find(ROASTcap3_left_tag);
                    
                    if(find_cp3 == string::npos) // to check if CAP3 contig belongs to ID (left or right CAP3)
                        find_cp3 = ID2.find(ROASTcap3_right_tag);
                    
                    if ((ID2 != ID) && (ID2.find(ROASTcap3_left_tag) != string::npos || ID2.find(ROASTcap3_right_tag) != string::npos) && (ID.find(ID2.substr(0, find_cp3)) != string::npos)) { //if cap3 found in ID, get cap3 derivative // if unique ID found in cap3 derivatives // it won't consider if contig ha split as multigene chimera
                        temp.push_back(ID2);
                        IDs2.erase(IDs2.begin() + ind2);
                        temp_ID.push_back(ID2);
                        //IDs2.remove(ID2);
                        ind2--;
                        it2--;
                    }
                    ind2++;
                }
                ind2 = 0;
                if (!(temp.empty())) {
                    //linked_contigs[ID].insert(temp);
                    linked_contigs.insert(make_pair(ID, temp));
                    temp.clear();
                    IDs.erase(IDs.begin() + ind);
                    temp_ID.push_back(ID);
                    //IDs.remove(ID);

                    ind--;
                    it--;
                }
            }
            ind++;
        }
        // cout << "made dictionary" << endl;
        bool left_done = false;
        bool right_done = false;
        int merge_count = 0;

        for (map<string, list<string> >::iterator it3 = linked_contigs.begin(); it3 != linked_contigs.end(); it3++) {
            string key = it3->first; // main ID

            final_fasta_seq = extract_fasta(key, in_assembly);
            // temp = it3->second;
            //list<string> cap3_assemblies(it3->second);
            for (list<string>::iterator it4 = it3->second.begin(); it4 != it3->second.end(); it4++) { //loop over cap3 derivatives

                string cap3_id = *it4;
                fasta_seq = extract_fasta(cap3_id, in_assembly);

                string cap3_ID_left =  ROASTcap3_left_tag;
                string cap3_ID_right =  ROASTcap3_right_tag;

                string cap3_ID_doubleleft = ROASTcap3_left_tag + "_" + ROASTcap3_left_tag;
                string cap3_ID_doubleright =  ROASTcap3_right_tag + "_" + ROASTcap3_right_tag;

                if (cap3_id.find(cap3_ID_left) != string::npos && left_done == false) { //cap3 left assembly?
                    final_fasta_seq = fasta_seq + string(Ns, 'N') + final_fasta_seq;
                    left_done = true;
                    merge_count++;

                } else if (cap3_id.find(cap3_ID_right) != string::npos && right_done == false) { //cap3 right assembly?
                    final_fasta_seq = final_fasta_seq + string(Ns, 'N') + fasta_seq;
                    right_done = true;
                    merge_count++;

                } else if (cap3_id.find(cap3_ID_doubleleft) != string::npos && left_done == true) { //cap3 left assembly?
                    final_fasta_seq = fasta_seq + string(Ns, 'N') + final_fasta_seq;
                    merge_count++;

                } else if (cap3_id.find(cap3_ID_doubleright) != string::npos && right_done == true) { //cap3 right assembly?
                    final_fasta_seq = final_fasta_seq + string(Ns, 'N') + fasta_seq;
                    merge_count++;
                }
            }

            if (ignore_short_seq == 0 || (ignore_short_seq > 0 && final_fasta_seq.size() > ignore_short_seq)) {

                if (change_header == true) {

                    ss << tool_name << "_" << contig_num;
                    string new_key = ss.str();
                    ss.str(string());
                    
                    out_file << ">" << new_key << endl << final_fasta_seq << endl;
                    contig_num++;
                    
                    size_t found = key.find("_and_");

                    if (found != string::npos) {

                        boost::replace_all(key, "_and_", ","); // replace all 'x' to 'y'
                        str_split(key, temp_ID2, ",");

                        summary << new_key << "\t";

                        for (int id = 0; id < temp_ID2.size(); id++) {

                            key = temp_ID2[id];
                            string key_temp;
                            int find = findNthOccur(key,'_', 1);

                            if (find > 0 && find <= key.size()) {

                                key_temp = key.substr(0, find);
                                summary << init_IDs_map.find(key_temp)->second << key.substr(find, key.size() - 1) << " , "; // to include original IDs plus modification tags
                                //cout << init_IDs_map.find(key_temp)->second << endl;
                               // cout << key.substr(find, key.size() - 1) << endl;

                            } else {
                                //cout << init_IDs_map.find(key)->second << " , ";
                                summary << init_IDs_map.find(key)->second << " , ";
                            }
                        }
                        summary << endl;
                        temp_ID2.clear();

                    } else {

                        string key_temp;
                        int find = findNthOccur(key, '_', 1);
                        
                         if (find > 0 && find <= key.size()) {
                            key_temp = key.substr(0, find);
                            summary << new_key << "\t" << init_IDs_map.find(key_temp)->second << key.substr(find, key.size() - 1) << endl; // to include original IDs plus modification tags

                        } else
                            summary << new_key << "\t" << init_IDs_map.find(key)->second << endl;
                    }

                } else { // don't change header name
                    out_file << ">" << key << endl << final_fasta_seq << endl;
                }
            } else if (ignore_short_seq > 0 && final_fasta_seq.size() <= ignore_short_seq) {

                del_seq << ">" << key << endl << final_fasta_seq << endl;
                summary  << "Removed: length <= " << ignore_short_seq << "\t" << ID << endl;
            }

            final_fasta_seq.empty();
            fasta_seq.empty();
            left_done = false;
            right_done = false;

        }
        in_file.open(in_assembly.c_str());
        //        cout << "fixed IDs" << endl;


        while (!getline(in_file, entry).eof()) {
            if (entry[0] == '>') {
                ID = entry.substr(1, entry.size() - 1);

            } else {
                if (find(temp_ID.begin(), temp_ID.end(), ID) == temp_ID.end()) { //current contig is not being merged

                    if (ignore_short_seq == 0 || (ignore_short_seq > 0 && entry.size() > ignore_short_seq)) {

                        if (change_header) {
                            ss << tool_name << "_" << contig_num;
                            string new_ID = ss.str();

                            ss.str(string());
                            out_file << ">" << new_ID << endl << entry << endl;
                            contig_num++;

                            size_t found = ID.find("_and_");
                            if (found != string::npos) {

                                boost::replace_all(ID, "_and_", ","); // replace all 'x' to 'y'
                                str_split(ID, temp_ID2, ",");

                                summary << new_ID << "\t";

                                for (int id = 0; id < temp_ID2.size(); id++) {

                                    ID = temp_ID2[id];
                                    string ID_temp;
                                    int find = findNthOccur(ID, '_', 1);

                                    if (find > 0 && find <= ID.size()) {

                                        ID_temp = ID.substr(0, find);
                                        summary << init_IDs_map.find(ID_temp)->second << ID.substr(find, ID.size() - 1) << " , "; // to include original IDs plus modification tags
                                       // cout << init_IDs_map.find(ID_temp)->second << ID.substr(find, ID.size() - 1) << " , ";

                                    } else
                                        summary << init_IDs_map.find(ID)->second << " , ";
                                   // cout << init_IDs_map.find(ID)->second << " , ";
                                }
                                summary << endl;
                                temp_ID2.clear();

                            } else {

                                string ID_temp;
                                int find = findNthOccur(ID, '_', 1);
                                
                                if (find > 0 && find <= ID.size()) {
                                    
                                    ID_temp = ID.substr(0, find);
                                    summary << new_ID << "\t" << init_IDs_map.find(ID_temp)->second << ID.substr(find, ID.size() - 1) << endl; // to include original IDs plus modification tags
                                   // cout << ID << "\t" << init_IDs_map.find(ID_temp)->second << ID.substr(find , ID.size() - 1) << endl;

                                } else
                                    summary << new_ID << "\t" << init_IDs_map.find(ID)->second << endl;
                            }
                        } else { // don't change header name
                            
                            out_file << ">" << ID << endl << entry << endl;
                        }
                    } else if (ignore_short_seq > 0 && entry.size() <= ignore_short_seq) {

                        del_seq << ">" << ID << endl << entry << endl;
                        summary << "Removed: length <= " << ignore_short_seq << "\t" << ID << endl;
                    }
                } else {
                    continue;
                }
            }
        }

        in_file.close();
        summary.close();
        del_seq.close();
        out_file.close();
        

        cout << "Number of cap3 asemblies merged to main contigs by Ns: " << merge_count << endl;
    }

  void update_contigIDs(string in_assembly, string out_assembly, string summary_file) {

        //avg_IS  = 120;
        std::size_t found = out_assembly.find_last_of("/\\");

        ofstream summary;
        summary.open(summary_file.c_str());
        
        ifstream in_file;
        in_file.open(in_assembly.c_str());

        string entry, ID, new_ID;
        int contig_num = 1;

        vector<string> IDs;
        stringstream ss;

        ofstream out_file;
        out_file.open(out_assembly.c_str());

        while (!getline(in_file, entry).eof()) {

            if (entry[0] == '>') {
                ID = entry.substr(1, entry.size() - 1);
            } else {
                ss << tool_name << "_" << contig_num;
                new_ID = ss.str();
                ss.str(string());
                summary << new_ID << "\t" << ID << endl;
                contig_num++;

                out_file << ">" << new_ID << endl << entry << endl;
            }
        }
        in_file.close();
        out_file.close();
        summary.close();
    }

    /*   void format_finalAssembly(string in_assembly, string out_assembly) {
           string entry, ID, ID2, fasta_seq, final_fasta_seq;
           map<string, list<string> >linked_contigs;
           list<string> temp;
           list <string> temp_ID;
           vector<string> IDs;
           ifstream in_file;
           int ind = 0, ind2 = 0;
           in_file.open(in_assembly.c_str());
           ofstream out_file;
           out_file.open(out_assembly.c_str());

           while (!getline(in_file, entry).eof()) { //save all IDs of fastafile
               if (entry[0] == '>') {
                   ID = entry.substr(1, entry.size() - 1);
                   //  cout << ID <<endl;
                   IDs.push_back(ID);
               }
           }
           in_file.close();
           cout << "got all IDs" << endl;
           vector<string> IDs2(IDs);
           for (std::vector<std::string>::iterator it = IDs.begin(); it != IDs.end(); ++it) {
               ID = *it;
               // cout << ID <<endl;
               if (ID.find("cap3") == string::npos) { // if cap3 not found , get unique ID
                   for (std::vector<std::string>::iterator it2 = IDs2.begin(); it2 != IDs2.end(); it2++) {
                       ID2 = *it2;
                       if ((ID2 != ID) && (ID2.find("cap3") != string::npos) && (ID2.find(ID) != string::npos)) { //if cap3 found in ID, get cap3 derivative // if unique ID found in cap3 derivatives
                           temp.push_back(ID2);
                           IDs2.erase(IDs2.begin() + ind2);
                           temp_ID.push_back(ID2);
                           //IDs2.remove(ID2);
                           ind2--;
                           it2--;
                       }
                       ind2++;
                   }
                   ind2 = 0;
                   if (!(temp.empty())) {
                       //linked_contigs[ID].insert(temp);
                       linked_contigs.insert(make_pair(ID, temp));
                       temp.clear();
                       IDs.erase(IDs.begin() + ind);
                       temp_ID.push_back(ID);
                       //IDs.remove(ID);

                       ind--;
                       it--;
                   }
               }
               ind++;
           }
           cout << "made dictionary" << endl;
           bool left_done = false;
           bool right_done = false;
        
           for (map<string, list<string> >::iterator it3 = linked_contigs.begin(); it3 != linked_contigs.end(); it3++) {
               string key = it3->first;

               final_fasta_seq = extract_fasta(key, in_assembly);
               // temp = it3->second;
               //list<string> cap3_assemblies(it3->second);
               for (list<string>::iterator it4 = it3->second.begin(); it4 != it3->second.end(); it4++) { //loop over cap3 derivatives
                
                   string cap3_id = *it4;
                   fasta_seq = extract_fasta(cap3_id, in_assembly);
                
                   string cap3_ID_left = key + "_cap3_left";
                   string cap3_ID_right = key + "_cap3_right";
                
                   if (cap3_id.find(cap3_ID_left) != string::npos && left_done == false) { //cap3 left assembly?
                       final_fasta_seq = fasta_seq + string(100, 'N') + final_fasta_seq;
                       left_done = true;
                    
                   } else if (cap3_id.find(cap3_ID_right) != string::npos && right_done == false) { //cape right assembly?
                       final_fasta_seq = final_fasta_seq + string(100, 'N') + fasta_seq;
                       right_done = true;

                   } else if (cap3_id.find("_cap3_left_cap3_left") != string::npos && left_done == true) { //cape left assembly?
                       final_fasta_seq = fasta_seq + string(100, 'N') + final_fasta_seq;

                   } else if (cap3_id.find("_cap3_right_cap3_right") != string::npos && right_done == true) { //cape right assembly?
                       final_fasta_seq = final_fasta_seq + string(100, 'N') + fasta_seq;
                   }
               }
               if (key.find("_and")) { // remove everything after _and in header
                   std::size_t found = key.find("_and");
                   out_file << ">" << key.substr(0, found) << endl << final_fasta_seq << endl;
               }
                   /*   else if(key.find("_cap3")) {
                          std::size_t found = key.find("_cap3");
                          out_file << ">" << key.substr(0, found) << endl << final_fasta_seq << endl;
                      }*/
    /*      else
              out_file << ">" << key << endl << final_fasta_seq << endl;

          final_fasta_seq.empty();
          fasta_seq.empty();
          left_done = false;
          right_done = false;

      }
      in_file.open(in_assembly.c_str());
      cout << "fixed IDs" << endl;


      while (!getline(in_file, entry).eof()) {
          if (entry[0] == '>') {
              ID = entry.substr(1, entry.size() - 1);
          } else {
              if (find(temp_ID.begin(), temp_ID.end(), ID) == temp_ID.end()) { //current contig is not being merged
                  if (ID.find("_and")) {
                      std::size_t found = ID.find("_and");
                      out_file << ">" << ID.substr(0, found) << endl << entry << endl;
                  } 

                  else
                      out_file << ">" << ID << endl << entry << endl;
              } else {
                  continue;
              }
          }
      }
      in_file.close();
      cout << "added remaining" << endl;
  } */
};
