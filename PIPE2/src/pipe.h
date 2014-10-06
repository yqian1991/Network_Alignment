/*
 * pipe.h
 *
 *  Created on: Sep 30, 2014
 *      Author: yqian33
 */

#ifndef PIPE_H_
#define PIPE_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <hash_map>

#define MATRIX_SIZE 200
//#define INVALID_SCORE -10
#define WINDOW_LEN 20
#define SPAM 35   //score of two substring, using Blosum 35
#define MATCH 10  //match number of substring 10

using namespace std;

namespace std{
   using namespace __gnu_cxx;
   struct str_hash{
	size_t operator()(const string& str) const{
		return  __stl_hash_string(str.c_str());

	}
   };
}

//double Blosum_matrix[MATRIX_SIZE][MATRIX_SIZE];//store the Blosum matrix value
/*********************************FILE CONSTANT*********************************************/
const string BLOSUM_FILE = "/home/yqian33/workspace/PIPE2/src/blosum80.txt";
const string PAM_FILE = "/home/yqian33/workspace/PIPE2/src/pam120.txt";
const string SEQ_FILE = "/home/yqian33/APPI/dataset/ProSeqBla/yeast/yeast559292.fasta";
const string PPI_FILE = "/home/yqian33/APPI/dataset/AlignMCLdataset/ppi/yeast.DIP.nif";
const string SCORE_FILE = "/home/yqian33/workspace/PIPE2/src/mat.txt";
/******************************************************************************************/

/******************************************************************************************/
const string pattern="\t";
vector<string> split(string str, string pattern);
void printMatrix(double blosum[][MATRIX_SIZE]);
/******************************************************************************************/


/******************************************SCORE FILE**************************************/
//Load score matrix
//Blosum file, not used
int readBlosum80(string infile);
double getBlosumScore(char a, char b);
double getBlosumScore(char a, char b, double m[][MATRIX_SIZE]);
//PAM120 file process, used
int readPam120(string infile);
double getPAMScore(char a, char b);
double getPAMScore(char a, char b, double m[][MATRIX_SIZE]);
/******************************************************************************************/


/*************************************DATA PRECOMPUTE**************************************/
//construct PPI network, save to hash_map
int interaction_Map(string infile, hash_map <string, vector<string>, str_hash > & inMap);//infile is interaction file
int check_inmap(string pname, map <string, vector<string> > & inMap);
int getSeqMap(string seqfile, hash_map<string, string, str_hash> &seqMap);
/******************************************************************************************/

/***********************************PREDICTION ENGINE**************************************/
//compute similarity of k-mer against a protein sequence using score matrix
int windows_similar(string ai, string pname);

//construct neighbor similar protein hash_map for every k-mer
int add_to_NbMap(hash_map<string, vector<string> , str_hash> &nbMap, string ai, vector<string> &prosecond);

//PPI predict engine
int pipe(string a, string b, hash_map<string, vector<string>, str_hash > &inMap, hash_map<string, string, str_hash> seqmap);//inMap is interaction_map
/******************************************************************************************/


#endif /* PIPE_H_ */
