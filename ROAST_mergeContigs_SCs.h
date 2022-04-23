/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ROAST_mergeContigs_SCs.h
 * Author: madiha
 *
 * Created on March 31, 2021, 9:23 PM
 */

#ifndef ROAST_MERGECONTIGS_SCS
#define ROAST_MERGECONTIGS_SCS

#include <string>


using namespace std;

class ROAST_mergeContigs_SCs {
public:
    ROAST_mergeContigs_SCs();
    ROAST_mergeContigs_SCs(const ROAST_mergeContigs_SCs& orig);
    virtual ~ROAST_mergeContigs_SCs();

    void FindFragmentsBySCs(string thread, string extended_assembly, string path_inter, string exe_path, int SC_cons_len_TH, int min_overlap_TH,  int blat_score_TH, int max_allowed_gaps);


private:

};

#endif /* ROAST_MERGECONTIGS_SCS_H */

