/*
 * AlignNemo.cpp
 *
 *  Created on: May 16, 2014
 *      Author: yqian33
 */
#include "anemo.h"
//#include "networkblast.h"

int mergeNif(string dir, string outfile)
{
		vector<string> files = vector<string>();
		getdir(dir,files);

		ofstream output;
		output.open(outfile.c_str());
		for(unsigned int i=2; i<files.size(); i++)
		{

			string fname = files[i].c_str();
			if(fname != "results_summary.txt" && fname != "alignment_graph.nif")
			{
				cout << fname << endl;
				ifstream input;
				string filename=dir+files[i];
				input.open(filename.c_str());

				if(!input.is_open())
				{
					cout << "file open error" << endl;
					return 0;
				}
				string solution="";
				string line;
				while( getline(input, line) )
				{
					vector<string> items = split(line, "\t");
					if(items.size()==4)
					{
						solution += items[0]+'\t'+items[1]+"\t";
					}
				}
				solution = solution.substr(0, solution.size()-1);
				output<<solution<<endl;
				input.close();
			}
		}
		output.close();
		cout<<"Merge nif files of AlignNemo"<<endl;
		return 1;
}


void testAnemo()
{
	string dir="/home/yqian33/APPI/software/AlignNemo-1.1/out_yeast-fly/";
	string outfile="/home/yqian33/APPI/software/AlignNemo-1.1/yeast-fly.solutions";
	mergeNif(dir, outfile);

}
