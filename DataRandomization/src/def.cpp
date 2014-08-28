/*
 * def.cpp
 *
 *  Created on: Aug 28, 2014
 *      Author: yqian33
 */
#include "def.h"
/*
 * infile: the input PPI interaction file
 * total: total interactions in infile
 * number: number of interactins you want to remove from infile
 * outfile: result interaction files.
 */
int remove_randomly(string infile, int total, int number, string outfile)
{
	ifstream input;
	input.open(infile.c_str());
	if( !input.is_open())
	{
		cout<<"open PPI interaction file"<<infile.c_str()<<" error\n";
		return 0;
	}

	ofstream output;
	output.open(outfile.c_str());

	int step = total/number;
	int start =random()%step;//the first line you want to remove, selected in random
	cout<<"start:"<<start<<" "<<"step:"<<step<<endl;
	int i=0;
	int j=0;//calculate the real number of output
	string line;
	while(getline(input, line)){
		if( (i-start)%step != 0)//not meet the requirement, then add to output
		{
			output<<line<<endl;
			j++;
		}
		i++;
	}
	cout<<"real number:"<<j<<endl;

	input.close();
	output.close();
	return 1;
}



