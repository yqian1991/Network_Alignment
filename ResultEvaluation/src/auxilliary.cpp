/*
 * auxilliary.cpp
 *
 *  Created on: May 6, 2014
 *      Author: yqian33
 */



#include "auxilliary.h"
#include "def.h"

using namespace std;

void statistic_known( )
{
	ifstream input;
	//input.open("/home/yqian33/APPI/software/AlignMCL-1.2/results/SS_scores/intra.ss.SimGIC.AlignMCL");
	input.open("/home/yqian33/APPI/software/AlignMCL-1.2/results/known_complexes_overlap/merged_PPIs/AlignMCL_overlap.txt");
	string line;
	int count = 0;
	int count1 = 0;
	int count2 = 0;
	while( getline(input, line))
	{
		vector<string> stritem = split(line, "\t");
		if (stritem[0] == "fly.all-yeast.all_fly")
		{
			count++;
			double f_score = atof(stritem[6].c_str());
			if (f_score - 0.3 >= 0 )
			{
				count1++;
			}
			if (f_score - 0.5 >= 0 )
			{
				count2++;
			}
		}

	}
	cout << "Matched solutions:"<< count<<endl;
	cout << "High quality solutions(F-index > 0.3):"<< count1<<endl;
	cout << "High quality solutions(F-index > 0.5):"<< count2<<endl;
}
