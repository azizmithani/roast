/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ROAST_extendContigs_SCs.h
 * Author: madiha
 *
 * Created on January 7, 2021, 12:24 PM
 */

#ifndef ROAST_EXTENDCONTIGS_SCS_H
#define ROAST_EXTENDCONTIGS_SCS_H
#include<map>
#include <string>
#include <vector>
#include <list>
#include <fstream> 
using namespace std;

class ROAST_extendContigs_SCs {
public:
    ROAST_extendContigs_SCs();
    ROAST_extendContigs_SCs(const ROAST_extendContigs_SCs& orig);
    virtual ~ROAST_extendContigs_SCs();
    void extend_bySoftclips(string fastafile, string ID_file, string BAM_FILE, string num, int SC_support_TH, string tool_name, int min_sc_reads, int sc_pos_from_corner, string exe_path);
    string consensus_seq_CAP3(std::vector<string>softclip_list, bool toleft, string thread, string exe_path, string PATH);
    string consensus_seq(std::vector<string>softclip_list, bool toleft);
    vector<string> extract_maxSC_block(vector<int> SC_base_pos, vector<string> SC_seq) ;
    int SC_start_pos(vector<int> SC_base_pos);
    
private:

};

#endif /* ROAST_EXTENDCONTIGS_SCS_H */

