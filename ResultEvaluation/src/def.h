/*
 * def.h
 *
 *  Created on: May 5, 2014
 *      Author: yqian33
 */

#ifndef DEF_H_
#define DEF_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <dirent.h>
#include <errno.h>

using namespace std;

const string pattern="\t";
vector<string> split(string str, string pattern);

int buildComplexMap(string infile, map <int, string> & proMap);
int buildSolMap(string infile, map <int, string> & proMap, int where);
void printMap(map <int, string> & proMap);

int overlap_Solution(map <int, string> & proMap, map <int, string> & complexMap, string outfile, string prefix, int complex_threshold);
int overlap_Complex(map <int, string> & proMap, map <int, string> & complexMap, string outfile, string prefix, int complex_threshold);
double F_index(string str1, string str2, int &count);
int find_overlap(string spe_pro, map <int, string> & proMap);
bool notExist(string single_pair, string spe_pro);

void clean_overlap(string infile);
int recovered(string str1, string str2);

int getdir (string dir, vector<string> &files);
#endif /* DEF_H_ */
