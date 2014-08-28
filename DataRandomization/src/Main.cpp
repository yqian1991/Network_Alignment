/*
 * Main.cpp
 *
 *  Created on: Aug 28, 2014
 *      Author: yqian33
 */

#include "def.h"

int main(void)
{
	string dir="/home/yqian33/workspace/PPIDataProcess/AddInteraction/src/yeast-test/";
	string infile=dir+"fly-yeast/yeast.DIP.nif";
	int total=22377;

	int number=10000;
	stringstream ss;
	ss<<number;

	int iter=10;
	int i=0;
	while(i<iter)
	{
		stringstream s1;
		s1<<i;
		string outfile=dir+"yeast_remove/"+"rm"+ss.str()+"_"+s1.str()+".nif";
		remove_randomly(infile, total, number, outfile);
		i++;
	}
}


