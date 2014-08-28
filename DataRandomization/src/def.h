/*
 * def.h
 *
 *  Created on: Aug 28, 2014
 *      Author: yqian33
 */

#ifndef DEF_H_
#define DEF_H_

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <algorithm>
#include <sstream>

using namespace std;

/*
 * remove interactions randomly from interaction file, eg.fly.DIP.nif
 */
int remove_randomly(string infile, int total, int number, string outfile);


#endif /* DEF_H_ */
