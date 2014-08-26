/*
 * networkblast.cpp
 *
 *  Created on: May 7, 2014
 *      Author: yqian33
 */

#include "networkblast.h"



string replaceChar(string str, char ch1, char ch2) {
  for (unsigned int i = 0; i < str.length(); ++i) {
    if (str[i] == ch1)
      str[i] = ch2;
  }

  return str;
}
/*
int getdir (string dir, vector<string> &files)
{
	 DIR* dirpath = opendir(dir.c_str());
	 dirent* pdir;
	 unsigned int num=0;
     while (pdir = readdir(dirpath)) {
    	 if(pdir->d_name[0]!='.'){
    		 num++;
    	 }
	 }
     //cout<<"file number:"<<num<<endl;
     string name;
     for(unsigned int i=2; i<=num+1; i++)
     {
    	 stringstream ss;
    	 ss<<i;
    	 name="";
    	 name = ss.str()+".sif";
    	 files.push_back(name);
    	 //cout<<name<<endl;

     }
     cout<<"file number:"<<files.size()<<endl;
     return 1;
}*/

int mergeSolutionsNB(string dir, string oufile)
{
	vector<string> files = vector<string>();
	getdir(dir,files);

	ofstream output;
	output.open(oufile.c_str());
	for(unsigned int i=0; i<files.size(); i++)
	{
		ifstream input;
		string filename=dir+files[i];
		input.open(filename.c_str());
		//cout << files[i].c_str()<< endl;
		if(!input.is_open())
		{
			cout << "file open error" << endl;
			return 0;
		}
		string solution="";
		string line;
		while( getline(input, line) )
		{
			vector<string> prostr = split(line, " ");
			//prostr[0], prostr[2]
			string pro1 = replaceChar(prostr[0], '|', '/');
			string pro2 = replaceChar(prostr[2], '|', '/');
			solution += pro1+"\t"+pro2+"\t";
		}
		solution = solution.substr(0, solution.size()-1);
		output<<solution<<endl;
		input.close();
	}
	output.close();
	cout<<"Merge solutions of NetworkBlast"<<endl;
	return 1;
}

string removeSpace(string str)
{
	int j=0;
	string result;
	for(unsigned int i=0; i<str.size(); i++)
	{
		 if(str[i]==' '){
		     ++j;
		 if(j==1){   // 这里改成双等号
		     //cout<<' ';
		     result +="\t";
		 }
		 if(i>1)
		     continue;
		 }else {
		     j=0;
		     result += str[i];
		 }
	}
	return result;
}
string map2UniprotAC(string name)
{
	string acname;
	locale settings;
    string uname;
    for(unsigned i = 0; i < name.size(); ++i)
    {
    		uname += (toupper(name[i], settings));
    }
    //cout << "uname:" <<uname << endl;

	ifstream input1;
	input1.open("/home/yqian33/APPI/ProSeqBla/idmapfile/YEAST_sysname_uniprotac_NEW.txt");
	string line;

	//cout << "name:" << name << endl;
	while( getline(input1, line) )
	{
			string line1=removeSpace(line);
			vector<string> item1 = split(line1, "\t");
			//cout <<"item1[0]:"<< item1[0]<<endl;
			if( uname==item1[0])
			{//remember item[1]: gene id
				//cout << "item1[1]:" << item1[1] << endl;
				acname=item1[1];
				return acname;
			}
			else
			{
				//acname="null";
				acname=uname;
			}

	}
	input1.close();
    //cout <<"END:"<< acname << endl;
    return acname;
}

int mipsComplexProcess(string infile,  map<string, string> & complex)
{
	ifstream input;
	input.open(infile.c_str());

	//ofstream output;
	//output.open(outfile.c_str());

	string line;
	while( getline(input, line) )
	{
		vector<string> items = split(line, "|");
		//just need items[0]: protein name and items[1]:catagory id
		vector<string> parts = split(items[1], ".");
		if( (items[1].find("550")==string::npos) && (parts.size()<4))
		{
			if ( existed(complex, items[1]))
			{
				//change the protein name
				//cout << "existed" <<endl;
				string str=complex[items[1]];
				//cout << "items[0]-pname:" << items[0] << endl;

				string ac = map2UniprotAC(items[0]);
				//cout << ac << endl;
				/*if( ac=="null")
				{
					ac=items[0];
				}*/
				//cout << ac << endl;
				str = str+" "+ac;
				//cout<< str << endl;
				complex[items[1]] = str;
			}
			else
			{
				//add
				//cout << "not existed" <<endl;
				string ac = map2UniprotAC(items[0]);
				/*if(ac=="null")
				{
					ac=items[0];
				}*/
				complex[items[1]]=ac;
			}
		}
	}

	input.close();
	cout << "size:"<<complex.size()<<endl;
	/*map<string, string>::iterator it = complex.begin();
	for(; it != complex.end(); it++)
	{
		output << it->first << "\t" << it->second << endl;
	}*/
	//output.close();
	return 1;
}

int existed(map<string, string> complex, string id)
{
	map<string, string>::iterator it = complex.begin();
	for(; it != complex.end(); it++)
	{
		//cout << it->first << endl;
		if (it->first == id )
		{
			//cout << it->second << endl;
			return 1;
		}
	}
	return 0;
}

int compare(string str1, string str2)
{


	int count=0;


	vector<string> solpros= split(str1, " ");
	vector<string> comPros= split(str2, " ");
	//cout <<"str11:"<< str1 << endl;
	//cout << str2 << endl;

	for (unsigned int i=0; i<solpros.size(); i++)
	{
		for (unsigned int j=0; j<comPros.size(); j++)
		{
			if( solpros[i] == comPros[j])
			{
				count++;
			}
		}
	}
	return count;
}
/*
 * a cluster is pure if it contained >=3 annotated proteins, and
 * at least half of them share the same annotation.
 */
int pureComplexVerify(string overlapFile, string complexFile, string outfile, int species)
{
	ifstream input;
	input.open(overlapFile.c_str());
	if (! input.is_open())
	{
		cout <<"open input file:" << overlapFile << "error!" << endl;
		return 0;
	}

	ofstream output;
	output.open(outfile.c_str());

	map <int, string> solMap = map <int, string> ();
	buildSolMap(overlapFile, solMap, species);

	map <string, string> complexMap = map <string, string> ();
	mipsComplexProcess(complexFile, complexMap);

	map<int, string>::iterator it = solMap.begin();

	int pureCount=0;
	bool sol_flag=false;
	int sol_count=0;

	for(; it != solMap.end(); it++)
	{
		vector<string> pro_item1 = split(it->second, " ");
		string str11="";
		for (unsigned int v=0; v<pro_item1.size()-1; v++)
		{
			if (notExist(pro_item1[v], str11))
			{
				//str11.append(pro_item1[v]+" ");
				str11 += pro_item1[v]+" ";
			}
		}

		//vector<string> solpros= split(it->second, " ");
		map<string, string>::iterator it1 = complexMap.begin();
		for(; it1 != complexMap.end(); it1++)
		{
			int count = compare(str11, it1->second);

			if (count>0)
			{
				output<<it->first+1<<"\t"<<it1->first<<"\t"<<count<<"\t"<<pro_item1.size()-1<<endl;
			}
			//cout << count << endl;
			if(count>=3)
			{
				pureCount++;
				//output<<it->first+1<<"\t"<<it1->first<<"\t"<<count<<endl;
				sol_flag=true;
			}
		}
		if(sol_flag)
		{
			sol_count++;
			sol_flag=false;
		}
	}
	double rate = (double) sol_count/solMap.size();
	cout<< rate << ":" << sol_count << " out of " << solMap.size() << endl;
	input.close();
	return 1;
}
