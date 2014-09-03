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

int addInteraction(int number, string infile1, string infile2, string outfile)
{
	ofstream ofs;
	ofs.open(outfile.c_str());

	ifstream ifs0;
	ifs0.open(infile1.c_str());

	ifstream ifs;
	ifs.open(infile2.c_str());


	string line;
	while(getline(ifs0, line))
	{
		ofs<<line<<endl;
	}
	ifs0.close();

	int i=0;
	cout<<"add more"<<endl;
	while(getline(ifs, line) && i<number)
	{
		vector<string> item = split(line, pattern);
		//ofs<<item[0]<<"\t"<<item[1]<<"\t"<<"1"<<endl;
		ofs<<line<<endl;
		i++;
	}

	ofs.close();
	ifs.close();

	return 1;



}

vector<string> split(string str, string pattern)
{
     string::size_type pos;
     vector<string> result;
     str+=pattern;//扩展字符串以方便操作
     int size=str.size();

     for(int i=0; i<size; i++)
     {
         pos=str.find(pattern,i);
         if(pos<size)
         {
             string s=str.substr(i,pos-i);
             result.push_back(s);
             i=pos+pattern.size()-1;
         }
     }
     return result;
}

