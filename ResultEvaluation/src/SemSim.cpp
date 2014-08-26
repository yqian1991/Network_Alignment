/*
 * semsim.cpp
 *
 *  Created on: May 13, 2014
 *      Author: yqian33
 */

#include "SemSim.h"
//#include "def.h"



/*
 * @params:solFile, the resulted alignment network file, each line in the file is a subnetwork solution
 * @params:goFile, Gene Association File, e.g."gene_association.goa_yeast"
 * @params:species, specify which species in the alignment you want to check
 * @params:outfile, output file contains the result, each line contains a protein name and all the GO terms annotated to it
 */
int compute_GOtermsPerProtein(string solFile, string gofile, string outfile, int species)
{
	ofstream output;
	output.open(outfile.c_str());

	map <int, string> solMap = map <int, string> ();
	buildSolMap(solFile, solMap, species);

	map<int, string>::iterator it = solMap.begin();
	for(; it != solMap.end(); it++)
	{

		vector<string> pro_item = split(it->second, " ");
		string str="";
		output<< "Solution " << it->first << endl;
		for (unsigned int v=0; v<pro_item.size()-1; v++)
		{
			if (notExist(pro_item[v], str))
			{
				string goterms;

				str += pro_item[v]+" ";

				//query GO terms for pro_item[v]
				ifstream input;
				input.open(gofile.c_str());

				string line;
				while ( getline(input, line))
				{
					vector<string> items = split(line, "\t");
					if (items[1] == pro_item[v])
					{
						//return GO terms id
						string gos = items[4];
						//cout << gos << endl;
						vector<string> go = split(gos, ":");
						goterms+=go[1]+" ";
					}
				}
				output<< pro_item[v] << "\t" << goterms <<endl;
				input.close();
			}
		}

	}
	output.close();
	cout << "Find Go terms for solutions done!" << endl;
	return 1;
}


/*
 * @params:solFile, the resulted alignment network file, each line in the file is a subnetwork solution
 * @params:goFile, go file here is the id mapping file from Uniprot, because it has Gene Association information. and easy for id mapping
 * @params:species, specify which species in the alignment you want to check
 * @params:outfile, output file contains the result, each line contains a protein name and all the GO terms annotated to it
 */
int compute_Goterms_MaWISh(string solFile, string gofile, string outfile, int species)
{
	ofstream output;
	output.open(outfile.c_str());

	cout <<"start..." << endl;

	map <int, string> solMap = map <int, string> ();
	buildSolMap(solFile, solMap, species);
	cout <<"solution map done, size="<< solMap.size() << endl;

	map<int, string>::iterator it = solMap.begin();
	for(; it != solMap.end(); it++)
	{
		vector<string> pro_item = split(it->second, " ");
		//cout << it->second << endl;
		string str="";

		output<< "Solution " << it->first << endl;
		for (unsigned int v=0; v<pro_item.size()-2; v++)
		{
			vector<string> gid = split(pro_item[v], "|");
			//cout << v << " "<< pro_item[v] <<" "<< gid[1] << endl;
			if (notExist(gid[1], str) && (gid[1] != "0"))
			{
				string goterms;
				str += gid[1]+" ";

				//query GO terms for gid[1]
				ifstream input;
				input.open(gofile.c_str());

				string line;
				while ( getline(input, line))
				{
					vector<string> items = split(line, "\t");
					//items[4] contains all the gi numbers
					if (items[4].find(gid[1]) != string::npos )//find GI
					{
						//cout<< "find" << endl;
						//return GO terms id
						string gos = items[6];
						//cout << gos << endl;
						//vector<string> go = split(gos, ":");
						goterms+=gos+" ";
						//cout<< gos << endl;
					}
				}
				output<< gid[1] << "\t" << goterms <<endl;
				input.close();
			}
		}
		cout<< "Solution " << it->first << "done" << endl;
	}
	output.close();
	cout << "Find Go terms for MaWISh solutions done!" << endl;
	return 1;
}

/*
 * @func:compute the intra species SS use SimGIC in fastSemSim
 * @params:solFile, the resulted alignment network file, each line in the file is a subnetwork solution
 * @params:dirname, specify a directory can store the temp file, each file hold one solution subnetwork
 * @params:species, specify the first or second species you want to operate
 */
int computeIntraSpeciesSimGIC(string solFile, int species, string dirname)
{
	map <int, string> solMap = map <int, string> ();
	buildSolMap(solFile, solMap, species);

	map<int, string>::iterator it = solMap.begin();
	int sol=0;
	for(; it != solMap.end(); it++)
	{
		stringstream ss;
		ss << sol;
		string filename = ss.str()+".txt";
		string filePath=dirname+filename;
		ofstream output;
		output.open(filePath.c_str());

		vector<string> pro_item = split(it->second, " ");
		string str="";

		for (unsigned int v=0; v<pro_item.size(); v++)
		{
			if (notExist(pro_item[v], str) )
			{
				str += pro_item[v]+" ";
				//cout<<pro_item[v]<<endl;
				output<<pro_item[v]<<endl;
			}
		}
		sol++;
		output.close();
	}

	vector<string> files = vector<string>();
	getdir(dirname,files);
	string outfile=solFile+".intra.ss";
	ofstream out;
	out.open(outfile.c_str());
	for(unsigned int i=2; i<files.size(); i++)
	{
		ifstream input;
		string filename=dirname+files[i];
		string outputfile=filename+".ss";
		string result=callFastSemSim(filename, outputfile, i-2, "list");
		out<<result<<endl;

	}
	out.close();
	return 1;
}

/*
 * @func: call fastSemSim actually
 * @params:filename, file contains the query list proteins, actually is the file in the directory in solution2queryList()
 * @params:number, file number
 * @params:qmode, lise or pair, parameters in fastSemSim
 */
string callFastSemSim(string filename,string output, int number, string qmode)
{
	string oboFile="/home/yqian33/APPI/software/AlignMCL-1.2/alignments/fastSemSim-0.8.1/test/data/Os/GeneOntology_filtered_2013.09.10.obo";

	string corpusFile1="/home/yqian33/APPI/software/AlignMCL-1.2/alignments/fastSemSim-0.8.1/test/data/ACs/gene_association.goa_fly_yeast";//for inter
	string corpusFile2="/home/yqian33/APPI/software/AlignMCL-1.2/alignments/fastSemSim-0.8.1/test/data/ACs/gene_association.goa_yeast";
	string corpusFile3="/home/yqian33/APPI/software/AlignMCL-1.2/alignments/fastSemSim-0.8.1/test/data/ACs/gene_association_bp.goa_fly";

	string progname="/home/yqian33/APPI/software/AlignMCL-1.2/alignments/fastSemSim-0.8.1/fastSemSim.sh ";
	string arg1="--o_type GeneOntology -o"+oboFile+" --o_file_format obo -a "+corpusFile2+" ";
	string arg2="--ac_type gaf2 --tss SimGIC --tmix BMA --query_ss_type obj --query_type SS ";
	//for BP:0008150   CC:0005575     MF:0003674
	string queryArgu="--query_input file --query_file  "+filename+" --query_mode "+qmode+"  --root GO:0008150 ";
	string outarg=" --output_file "+output;

	string command=progname+arg1+arg2+queryArgu+outarg;

	system(command.c_str());
	cout<<"\n-----------------------------------------------\n";

	ifstream input;
	input.open(output.c_str());

	string line;
	int pair_count=0;
	int mean_count=0;
	double sum_ss=0.0;
	while ( getline(input, line) )
	{

		vector<string> items = split(line, "\t");
		if( items[0]!= items[1])
		{
			//cout<<line<<endl;
			pair_count++;
			//if( items[2] != "None" && items[2] != "0.0")
			if( items[2] != "None")
			{
				mean_count++;
				sum_ss += atof(items[2].c_str());
			}
		}
	}
	cout<<"-----------------------------------------------\n";
	double ss = sum_ss/pair_count;
	cout<<number<<"\t"<<ss<<"\t"<<mean_count<< "\t" <<pair_count<<endl;

	ostringstream sstream;
	sstream<<number<<"\t"<<ss<<"\t"<<mean_count<< "\t" <<pair_count;
	string res=sstream.str();
	input.close();
	return res;
}

/*
 * @func:compute the inter species SS using fastSemSim SimGIC
 * @params:solFile, the resulted alignment network file, each line in the file is a subnetwork solution
 * @params:dirname, the directory contains all the solutions, one file holds one solution
 */
int computeInterSpeciesSimGIC(string solFile, string dirname)
{
	map <int, string> solMap1 = map <int, string> ();
	buildSolMap(solFile, solMap1, 0);

	map <int, string> solMap2 = map <int, string> ();
	buildSolMap(solFile, solMap2, 1);

	cout<<"Solution Map For Two Species Constructed."<<endl;

	map<int, string>::iterator it1 = solMap1.begin();
	map<int, string>::iterator it2 = solMap2.begin();
	int sol=0;
	for(; it1 != solMap1.end(), it2 != solMap2.end(); it1++, it2++)
	{
		stringstream ss;
		ss << sol;
		string filename = ss.str()+".txt";
		string filePath=dirname+filename;
		ofstream output;
		output.open(filePath.c_str());

		vector<string> pro_item1 = split(it1->second, " ");
		string str1="";//used for clean repeat proteins
		for (unsigned int v=0; v<pro_item1.size(); v++)
		{
			if (notExist(pro_item1[v], str1) )
			{
				str1 += pro_item1[v]+" ";
			}
		}

		vector<string> pro_item2 = split(it2->second, " ");
		string str2="";//used for clean repeat proteins
		for (unsigned int v=0; v<pro_item2.size(); v++)
		{
			if (notExist(pro_item2[v], str2) )
			{
				str2 += pro_item2[v]+" ";
			}
		}
		//str1, str2 to file
		vector<string> p1 = split(str1, " ");
		vector<string> p2 = split(str2, " ");
		for(unsigned int i=0; i<p1.size()-1; i++)
		{
			for(unsigned int j=0; j<p2.size()-1; j++)
			{
				output<<p1[i]<<"\t"<<p2[j]<<"\n";//used for qmode pair in fastSemSim
			}
		}
		sol++;
		output.close();
	}

	vector<string> files = vector<string>();
	getdir(dirname,files);
	string outfile=solFile+".inter.ss";
	ofstream out;
	out.open(outfile.c_str());
	for(unsigned int i=2; i<files.size(); i++)
	{
		string filename=dirname+files[i];
		string outputfile=filename+".ss";
		string result=callFastSemSim(filename, outputfile, i-2, "pairs");
		out<<result<<endl;

	}
	out.close();

	return 1;
}

/*
 * Extract BP terms of a GOA file
 */
int pickBP(string gotermsFile, string bpFile, string out)
{
	ifstream input1;
	input1.open(bpFile.c_str());

	vector<string> bp;
	vector<string>::iterator iter;
	string line;
	while( getline(input1, line))
	{
		vector<string> item = split(line, " ");
		//if (item[0] =="Biological Process" || item[0] =="Molecular Function" )
		iter=bp.end();
		//vector<string> gos = split(item[1], ":");
		bp.insert(iter, item[0]);
	}
	cout<<"Read and Store Biological Process GO terms..."<<endl;
	input1.close();

	ifstream input2;
	input2.open(gotermsFile.c_str());

	ofstream output;
	output.open(out.c_str());

	int count=0;
	int number=0;
	while( getline(input2, line) )
	{
		if(count==0)
		{
			output<<line<<"\n";
		}
		else if(count>0)
		{
			vector<string> items = split(line, "\t");
			string go = items[4];
			if( find(bp.begin(), bp.end(), go) != bp.end() )
			{
				number++;
				//cout<<"find"<<endl;
				output<<line<<"\n";
			}
		}
		count++;
	}
	cout<<number<<" terms belong to BP"<<endl;
	input2.close();
	output.close();
	return 1;
}

/*
 * Extract all the BP terms assigned to each solution
 */
int pickBPtermsOfSolutions(string gotermsFile, string bpFile,  string out)
{
	ifstream input1;
	input1.open(bpFile.c_str());

	vector<string> bp;
	vector<string>::iterator iter;
	string line;
	while( getline(input1, line))
	{
		vector<string> item = split(line, " ");
		iter=bp.end();
		bp.insert(iter, item[0]);
	}
	cout<<"Reading and Store Biological Process GO terms..."<<endl;
	input1.close();

	ifstream input2;
	input2.open(gotermsFile.c_str());

	ofstream output;
	output.open(out.c_str());

	while( getline(input2, line))
	{
		//cout<<line<<endl;

		vector<string> items= split(line, "\t");
		if(items[0].find("Solution")!=string::npos)
		{
			output<<line<<"\n";
		}
		else
		{
			output<<items[0]<<"\t";

			string termStr="";
			vector<string> terms = split(items[1], " ");

			for(unsigned int i=0; i<terms.size()-1; i++)
			{
				//cout<<terms[i]<<endl;
				if( find(bp.begin(), bp.end(), "GO:"+terms[i])!= bp.end())
				{
					termStr=termStr +terms[i]+" ";
					//cout<<termStr<<endl;
					output<<terms[i].c_str()<<" ";
				}
			}
			output<<"\n";
		}
	}
	input2.close();
	output.close();
	return 1;

}
void test()
{
	string bpFile="/home/yqian33/APPI/software/AlignMCL-1.2/alignments/fastSemSim-0.8.1/test/data/ACs/GO_BF.txt";


	/*string gotermsFile="/home/yqian33/APPI/software/AlignMCL-1.2/alignments/fastSemSim-0.8.1/test/data/ACs/gene_association.goa_fly";
	string out="/home/yqian33/APPI/software/AlignMCL-1.2/alignments/fastSemSim-0.8.1/test/data/ACs/gene_association_bp.goa_fly";
	if( pickBP(gotermsFile, bpFile, out))
	{
		cout<<"done"<<endl;
	}*/

	string gotermsFile="/home/yqian33/APPI/software/MAWISH/data/alignment/fly-yeast-i2d/fly-yeast_yeast-goterms.txt";
	string out        ="/home/yqian33/APPI/software/MAWISH/data/alignment/fly-yeast-i2d/fly-yeast_yeast-goterms-bp.txt";

	if( pickBPtermsOfSolutions(gotermsFile, bpFile, out))
	{
		cout<<"done"<<endl;
	}
}
