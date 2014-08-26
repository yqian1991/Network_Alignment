/*
 * networkblast.h
 *
 *  Created on: May 7, 2014
 *      Author: yqian33
 */

#ifndef NETWORKBLAST_H_
#define NETWORKBLAST_H_

#include <iostream>
#include <cstdlib>
#include <string>
#include <cstring>
#include <vector>
#include <algorithm>
#include <locale>
#include "auxilliary.h"
#include "def.h"


using namespace std;



string replaceChar(string str, char ch1, char ch2);
//int getdir (string dir, vector<string> &files);
int mergeSolutionsNB(string dir, string outfile);

int pureComplexVerify(string overlapFile, string complexFile, string outfile, int species);

string removeSpace(string str);
string map2UniprotAC(string name);

int mipsComplexProcess(string infile, map<string, string> & complex);
int existed(map<string, string> complex, string id);



#endif /* NETWORKBLAST_H_ */
