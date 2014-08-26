/*
 * netaligner.cpp
 *
 *  Created on: Jun 18, 2014
 *      Author: yqian33
 */

#include "netaligner.h"


int mergeSolutions(string dir, string oufile)
{
	vector<string> files = vector<string>();
	getdir(dir,files);
	sort(files.begin(), files.end());

	ofstream output;
	output.open(oufile.c_str());

	for(unsigned int i=2; i<files.size(); i++)
	{
		ifstream input;

		string filename = files[i];
		if(filename[0]=='a')
		{
			string filename=dir+files[i];
			input.open(filename.c_str());
			//cout << files[i]<< endl;
			if(!input.is_open())
			{
				cout << "file open error" << endl;
				return 0;
			}

			string line="";
			int count=0;
			string first="";
			vector<string>oneSolution;
			vector<string>::iterator iter;
			while( getline(input, line) )
			{
				if( count>0 )
				{
					vector<string>items = split(line, "\t");
					iter = oneSolution.end();
					oneSolution.insert(iter, items[0]);
					iter = oneSolution.end();
					oneSolution.insert(iter, items[1]);

				}
				count++;
			}
			sort( oneSolution.begin(), oneSolution.end() );
			iter = unique( oneSolution.begin(), oneSolution.end() );
			oneSolution.erase(iter, oneSolution.end());
			for(unsigned int i =0; i<oneSolution.size(); i++)
			{
				output<<oneSolution[i]<<"\t";
			}
			output<<"\n";
			input.close();
		}

	}
	output.close();
	cout<<"Merge solutions of NetworkBlast"<<endl;
	return 1;
}
