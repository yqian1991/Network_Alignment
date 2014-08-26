/*
 * SemSim.h
 *
 *  Created on: May 13, 2014
 *      Author: yqian33
 */

#ifndef SEMSIM_H_
#define SEMSIM_H_

#include <iostream>
#include <string>
#include <cstdlib>
//#include "auxilliary.h"
#include "networkblast.h"


using namespace std;

//gotermsFile is the solution subnetworks with GO terms annotated to each protein in the subnetwork
//bpFile is the GO terms which specified their categories(BP, CC, MF)
//1. pick out BP terms and annotate to solutions
//2. pick out BP terms and calcualte Semantic Similarity
int pickBP(string gotermsFile, string bpFile,  string out);

//pick from annotated solution subnetwork file
int pickBPtermsOfSolutions(string gotermsFile, string bpFile,  string out);

int compute_GOtermsPerProtein(string solFile, string gofile, string outfile, int species);

int compute_Goterms_MaWISh(string solFile, string gofile, string outfile, int species);

//for SS calculation
int computeInterSpeciesSimGIC(string solFile, string dirname);

int computeIntraSpeciesSimGIC(string solFile, int species, string dirname);

string callFastSemSim(string filename,string output, int number, string qmode);



void test();

#endif /* SEMSIM_H_ */
