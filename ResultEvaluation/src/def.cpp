/*
 * def.cpp
 *
 *  Created on: May 5, 2014
 *      Author: yqian33
 */
#include "def.h"

/*
 * @func:store the known complexes to map
 * @params:infile, the known complexes file
 * @params:proMap, the map you want to store complexes into
 */
int buildComplexMap(string infile, map <int, string> & proMap)
{
	ifstream input;
	input.open(infile.c_str());
	if (!input.is_open())
	{
		cout <<"Open File"<<infile<<" Error"<<endl;
		return 0;
	}
	string line;
	unsigned int pid = 0;
	string protein="start";
	while ( getline(input, line))
	{

		vector<string> strsplit = split(line, pattern);
		unsigned int cid = atoi(strsplit[0].c_str());

		string pro = strsplit[1];
		//cout << "pro: "<<pro << endl;
		if( pid != cid)//write a new protein
		{
			//cout << protein << endl;
			proMap[pid]=protein;
			//protein = " ";
			//protein.append("\t"+pro);
			protein=pro;
		}
		else
		{
			protein.append(" "+pro);
		}
		pid = cid;
	}
	input.close();

	return 1;

}
/*
 * @func:store the resulted alignment networks to map
 * @params:infile, the result alignment networks file
 * @params:proMap, the map you want to store networks into
 * @params:where, specify the first or the second species in the alignment network you want to extract
 */
int buildSolMap(string infile, map <int, string> & proMap, int where)
{
	ifstream input;
	input.open(infile.c_str());
	if (!input.is_open())
	{
		cout <<"Open File"<<infile<<" Error"<<endl;
		return 0;
	}
	string line;
	int sol_num=0;
	string protein="";
	//cout<<"here"<<endl;
	while ( getline(input, line))
	{
		//cout<< line << endl;
		vector<string> strsplit = split(line, pattern);

		for(unsigned int i=0; i<strsplit.size(); i++)//as some solution files（MaWISh） contain blank, so i should up to strsplit.size()-1
		{
			//if (strsplit[i].find("/") != std::string::npos)//if so
			//{
				vector<string> proName = split(strsplit[i], "/");
				//cout<< proName[where] << endl;
				protein += proName[where]+" ";
			//}
		}

		proMap[sol_num++]=protein;
		//cout << sol_num-1 << " "<< protein << endl;
		protein="";

	}
	input.close();

	return 1;
}
/*
 * consistent with the result in known_complexes_overlap/single_PPIs/AlignMCL_overlap.txt
 * @params:proMap, is the map for solutions resulted from network alignment, one key map to one subnetwork
 * @params:complexMap, is the map for the known complexes, one key map to one complex
 * @params:outfile, output contains information about overlap
 * @params:prefix, prefix at the begining of each line in output file
 * @params:complex_threshold, specify the smallese number of size of which complexes you will check
 */

int overlap_Solution(map <int, string> & proMap, map <int, string> & complexMap, string outfile, string prefix, int complex_threshold)
{
	cout<<"Run overlap from Solutions..."<<endl;
	ofstream output;
	output.open(outfile.c_str());

	int total=0;
	int highsol_count=0;
	int highersol_count=0;

	map<int, string>::iterator pro_it = proMap.begin();
	for(; pro_it != proMap.end(); pro_it++)
	{
		int match_gt_2=0;
		int complexId, complexNum;
		int proId = pro_it ->first;
		vector<string> pro_item0 = split(pro_it->second, " ");
		string str2="";
		for (unsigned int v=0; v<pro_item0.size()-1; v++)
		{
			if (notExist(pro_item0[v], str2))
			{
				str2.append(pro_item0[v]+" ");
			}
		}
		vector<string> pro_item = split(str2, " ");
		int proNum = pro_item.size()-1;

		double f_score=0.00;
		int tp_count=0;

		//match_gt_2=0;
		map<int, string>::iterator comp_it = complexMap.begin();
		for(; comp_it != complexMap.end(); comp_it++)
		{
			vector<string> compitem = split(comp_it->second, " ");
			complexNum = compitem.size()-1;
			if(complexNum>=complex_threshold)
			{
				int tp_num=0;
				double rawf = F_index(comp_it->second, str2, tp_num);
				if( rawf >0)
				{
					match_gt_2++;
					if (rawf-f_score>0)
					{
						f_score = rawf;
					}
					complexId = comp_it->first;

					tp_count=tp_num;
					output << prefix << "\t"<< proId << "\t" << complexId << "\t" << proNum << "\t" <<  complexNum << "\t" << tp_count << "\t" << rawf << endl;

				}
			}
		}
		if(match_gt_2>1)
		{
			total++;
		}
		if( f_score-0.3>=0)
		{
			highsol_count++;
		}
		if( f_score-0.5>=0)
		{
			highersol_count++;
		}
	}

	cout << "Matched solutions:"<< total<<endl;
	cout << "High quality solutions(F-index >= 0.3):"<< highsol_count<<endl;
	cout << "High quality solutions(F-index >= 0.5):"<< highersol_count<<endl;
	output << "Matched solutions:"<< total<<endl;
	output << "High quality solutions(F-index >= 0.3):"<< highsol_count<<endl;
	output << "High quality solutions(F-index >= 0.5):"<< highersol_count<<endl;
	output.close();
	return 1;
}

/*
 * consistent with the result in known_complexes_overlap/single_PPIs/AlignMCL_overlap.txt
 * infile is known complexes file.
 * consistent with the result in known_complexes_overlap/single_PPIs/AlignMCL_overlap.txt
 * @params:proMap, is the map for solutions resulted from network alignment, one key map to one subnetwork
 * @params:complexMap, is the map for the known complexes, one key map to one complex
 * @params:outfile, output contains information about overlap
 * @params:prefix, prefix at the begining of each line in output file
 * @params:complex_threshold, specify the smallese number of size of which complexes you will check
 */
int overlap_Complex(map <int, string> & proMap, map <int, string> & complexMap, string outfile, string prefix, int complex_threshold)
{
	cout<<"Run overlap from complexes..."<<endl;
	ofstream output;
	output.open(outfile.c_str());

	//for every complex , search matched sols
	map<int, string>::iterator comp_it = complexMap.begin();
	comp_it++;

	int highsol_count=0;
	int highersol_count=0;
	int scr=0;
	for(; comp_it != complexMap.end(); comp_it++)
	{

		//cout << comp_it->second<< endl;
		int proId, proNum;
		int complexId = comp_it ->first;
		vector<string> compitem = split(comp_it->second, " ");
		int complexNum = compitem.size();

		if(complexNum>=complex_threshold)
		{
			double f_score=0.00;
			int tp_count=0;

			//proteins in comp_it->second;
			map<int, string>::iterator pro_it = proMap.begin();
			for(; pro_it != proMap.end(); pro_it++)
			{
				//cout << pro_it->second<< endl;
				vector<string> pro_item = split(pro_it->second, " ");
				string str2="";
				for (unsigned int v=0; v<pro_item.size()-1; v++)
				{
					if (notExist(pro_item[v], str2))
					{
						str2.append(pro_item[v]+" ");
					}
				}
				int tp_num=0;

				vector<string> str2item=split(str2, " ");
				if (str2item.size()-1>1)
				{
					double rawf = F_index(comp_it->second, str2, tp_num);
					if( rawf-f_score >0)
					{
						f_score = rawf;
						proId = pro_it->first;

						vector<string> proitem = split(str2, " ");
						proNum = proitem.size()-1;
						tp_count=tp_num;
					}
				}
			}
			//
			output << prefix << "\t"<< complexId << "\t" << proId << "\t" << complexNum << "\t" <<  proNum << "\t" << tp_count << "\t" << f_score << endl;
			if( f_score-0.3>=0)
			{
				highsol_count++;
			}
			if( f_score-0.5>=0)
			{
				highersol_count++;
			}
		}
		else if( complexNum<complex_threshold && complexNum>1 )
		{//calculate the small recovered complexes
			bool scrFlag=false;
			map<int, string>::iterator pro_it_src = proMap.begin();
			for(; pro_it_src != proMap.end(); pro_it_src++)
			{
				vector<string> proteinItem = split(pro_it_src->second, " ");
				string str3="";
				for (unsigned int v=0; v<proteinItem.size()-1; v++)
				{
					if (notExist(proteinItem[v], str3))
					{
						str3.append(proteinItem[v]+" ");
					}
				}

				vector<string> protein_rm = split(str3, " ");
				proNum =  protein_rm.size()-1;
				if ( proNum < 20)
				{
					if ( recovered(comp_it->second, str3) )
					{
						scrFlag=true;
					}
				}

			}
			if (scrFlag)
			{
				scr++;
			}
		}
	}

	//cout << "Matched solutions:"<< sol_count<<endl;
	cout << "High quality solutions(F-index > 0.3):"<< highsol_count<<endl;
	cout << "High quality solutions(F-index > 0.5):"<< highersol_count<<endl;
	cout << "Recovered Small Complexes:"<< scr<<endl;

	output << "High quality solutions(F-index > 0.3):"<< highsol_count<<endl;
	output << "High quality solutions(F-index > 0.5):"<< highersol_count<<endl;
	output << "Recovered Small Complexes:"<< scr<<endl;
	output.close();
	return 1;
}
/*
 * @func: at least two items in str1 should exist in str2 return 1, else return 0
 */
int recovered(string str1, string str2)
{
	int count=0;
	int flag=0;

	vector<string> str1_item = split(str1, " ");
	vector<string> str2_item = split(str2, " ");
	for(unsigned int i=0; i<str1_item.size(); i++)
	{
		for(unsigned int j=0; j<str2_item.size(); j++)
		{
				if( str1_item[i]== str2_item[j])
				{
					count++;
				}
		}
	}
	if(count>1)
	{
		flag = 1;
	}
	else{
		flag = 0;
	}
	return flag;
}
/*
 * @func:compute the F_index(recall and presicion of overlapping)
 */
double F_index(string str1, string str2, int &count)
{
	double fscore=0.00;
	vector<string> str1_item = split(str1, " ");
	vector<string> str2_item = split(str2, " ");
	//int count=0;
	for(unsigned int i=0; i<str2_item.size()-1; i++)
	{
		for(unsigned int j=0; j<str1_item.size(); j++)
		{
			if (str2_item[i]== str1_item[j])
			{
				count++;
			}
		}
	}

	//if (count != 0)
	//{
		double pi = (double)count/(str2_item.size()-1);
		double ro = (double)count/str1_item.size();
		fscore = (double)2*pi*ro/(pi+ro);
	//}
	return fscore;
}

/*
 * @func:check whether single_pair is in spe_pro or not
 * @params:single_pair, the protein str you want to check
 * @params:spe_pro, the repeat masked proteins string, you want to query against, and add
 */
bool notExist(string single_pair, string spe_pro)
{
	bool flag = true;
    vector<string> str = split(spe_pro, " ");
    for (unsigned int i=0; i<str.size(); i++)
    {
    	if(single_pair == str[i])
    	{
    		flag = false;
    	}
    }
	return flag;

}

int find_overlap(string spe_pro, map <int, string> & proMap)
{
	int count=0;

	map<int, string>::iterator it = proMap.begin();
	it++;
	for(; it != proMap.end(); it++)
	{
		vector<string> subject = split(it->second, " ");
		vector<string> query = split(spe_pro, " ");
		for(unsigned int i=0; i<query.size(); i++)
		{
			for(unsigned int j=0; j<subject.size(); j++)
			{
				if (query[i] == subject [j])
				{
					count++;

				}
			}
		}

	}
	return count;

}

void clean_overlap(string infile)
{
	ifstream input;
	input.open(infile.c_str());

	string line;
	while( getline(input, line))
	{
		vector<string> str = split(line, " ");
		//str[1] as key, line as value
		//string value = str[0]+
	}
}
void printMap(map <int, string> & proMap)
{
	map<int, string>::iterator it = proMap.begin();
	for(; it != proMap.end(); it++)
	{
		cout<<it->first<<'\t'<<it->second<<endl;
	}
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
}
