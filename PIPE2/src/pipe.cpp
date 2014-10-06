/*
 * pipe.cpp
 *
 *  Created on: Sep 30, 2014
 *      Author: yqian33
 */
#include "pipe.h"

double Blosum_matrix[MATRIX_SIZE][MATRIX_SIZE];//store the Blosum matrix value
double pam_matrix[MATRIX_SIZE][MATRIX_SIZE];//store the PAM120 matrix value

vector<string> split(string str, string pattern)
{
     string::size_type pos;
     vector<string> result;
     str+=pattern;//扩展字符串以方便操作
     unsigned int size=str.size();

     for(unsigned int i=0; i<size; i++)
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

void printMatrix(double blosum[][MATRIX_SIZE])
{
	for(int i=0; i<MATRIX_SIZE; i++)
	{
		for(int j=0; j<MATRIX_SIZE; j++)
		{
			cout<<*(*(blosum+i)+j)<<" ";
		}
		cout<<endl;
	}
}
/*
 * @name: readBlosum80
 * @return: 2D array that store the score of Blosum80
 */
int readBlosum80(string infile)//, double matrix[][MATRIX_SIZE])
{
	ifstream input;
	input.open(infile.c_str());
	if ( !input.is_open())
	{
		cout<<infile<<" open error"<<endl;
		return 0;
	}
	string line;

	getline(input, line);
	vector<string> acid_name = split(line.substr(1), "  ");
	//cout<<acid_name.size()<<endl;
	int acid_index[25]={0};
	for (unsigned int i=0; i<acid_name.size(); i++)
	{
		const char *ch = acid_name[i].c_str();
		acid_index[i] = int(*ch);
		//cout<<*ch<<" "<<acid_index[i]<<endl;
	}
	int i=1;
	int index=1;
	while( getline(input, line) )
	{
		double score=0.0;
		//process the line to get score
		vector<string> score_str = split(line.substr(1), " ");
		index=1;
		for(unsigned int j=0; j<score_str.size(); j++)
		{
			if (score_str[j] == " " or score_str[j] == "  " or score_str[j] == "" )
				continue;
			score = atof(score_str[j].c_str());
			//cout<<i<<" "<<index<<" "<<acid_index[i]<<" "<<acid_index[index]<<" "<<score<<endl;
			Blosum_matrix[acid_index[i]][acid_index[index++]] = score;
		}
		i++;
	}
	input.close();
	return 1;
}
double getBlosumScore(char a, char b)
{
	return Blosum_matrix[int(a)][int(b)];
}
double getBlosumScore(char a, char b, double m[][MATRIX_SIZE])
{

	return m[int(a)][int(b)];
}

/*
 * @FUNCTION: Load the score matrix of PAM120
 */
int readPam120(string infile)//, double matrix[][MATRIX_SIZE])
{
	ifstream input;
	input.open(infile.c_str());
	if ( !input.is_open())
	{
		cout<<infile<<" open error"<<endl;
		return 0;
	}
	string line;

	getline(input, line);
	vector<string> acid_name = split(line.substr(1), "  ");
	//cout<<acid_name.size()<<endl;
	int acid_index[25]={0};
	for (unsigned int i=0; i<acid_name.size(); i++)
	{
		const char *ch = acid_name[i].c_str();
		acid_index[i] = int(*ch);
		//cout<<*ch<<" "<<acid_index[i]<<endl;
	}
	int i=1;
	int index=1;
	while( getline(input, line) )
	{
		double score=0.0;
		//process the line to get score
		vector<string> score_str = split(line.substr(1), " ");
		index=1;
		for(unsigned int j=0; j<score_str.size(); j++)
		{
			if (score_str[j] == " " or score_str[j] == "  " or score_str[j] == "" )
				continue;
			score = atof(score_str[j].c_str());
			//cout<<i<<" "<<index<<" "<<acid_index[i]<<" "<<acid_index[index]<<" "<<score<<endl;
			pam_matrix[acid_index[i]][acid_index[index++]] = score;
		}
		i++;
	}
	input.close();
	return 1;
}

/*
 * @FUNCTION: Get the score of two amino acid
 */
double getPAMScore(char a, char b)
{
	return pam_matrix[int(a)][int(b)];
}
double getPAMScore(char a, char b, double m[][MATRIX_SIZE])
{
	return m[int(a)][int(b)];
}

/*
 * @interaction__Map
 * @return: inMap is interaction map(namely protein network), key is protein name, value is its neighbors
 */
int interaction_Map(string infile, hash_map <string, vector<string>, str_hash > & inMap)
{
	ifstream input;
	input.open(infile.c_str());
	if ( !input.is_open())
	{
		cout<<infile<<" open error"<<endl;
		return 0;
	}
	string line="";
	while( getline(input, line))
	{
		vector<string> pname= split(line, "\t");
		//add protein 1
		if (inMap.find(pname[0]) != inMap.end())
		{
			inMap[pname[0]].push_back(pname[1]);
		}
		else
		{
			vector<string> nb;
			nb.push_back(pname[1]);
			inMap[pname[0]] = nb;
		}
		//add protein 2
		if (inMap.find(pname[1]) != inMap.end())
		{
			inMap[pname[1]].push_back(pname[0]);
		}
		else
		{
			vector<string> nb;
			nb.push_back(pname[0]);
			inMap[pname[1]].push_back(pname[0]);
		}

	}
	cout<<"Interaction Graph size: "<<inMap.size()<<endl;
	input.close();
	return 1;
}

int check_inmap(string pname, hash_map <string, vector<string>, str_hash > & inMap)
{
	return inMap[pname].size();
}
/*
 * @Input: interaction graph, inMap
 * @Input: a , b are two protein sequences
 * @Output: matrix |A|*|B|
 * @Estimate: Need 50min to judge two proteins.
 */
int pipe(string a, string b, hash_map<string, vector<string>, str_hash > &inMap, hash_map<string, string, str_hash> seqmap)
{
	string pname="";//used in the iteration ,store protein names in inMap
	string seq="";//used in the iteration ,store sequences of pname

	int aLen = a.length();
	int bLen = b.length();
	string ai="";//store the window sequence of ai
	string bi="";//store the window sequence of bi

	//prepare sequence map key is pname, value is sequence
	//map<string, string> seqmap;
	//getSeqMap(SEQ_FILE, seqmap);//this can do out of the function

	hash_map<string, vector<string>, str_hash > nbMap;//find similar of ai, and add neighbors here
	hash_map<string, vector<string>, str_hash >::iterator nbMap_it;//iterator of nbMap
	hash_map<string, vector<string>, str_hash >::iterator pro_it;//iterator of interaction map

	for(int i=0; i<=aLen-WINDOW_LEN; i++){
		if (i==0){
			ai = a.substr(i, WINDOW_LEN);
		}
		else{
			ai = ai.erase(0, 1)+a[i+WINDOW_LEN-1];
			//cout<<ai<<endl;
		}
		//cout<<"ai:"<<ai<<endl;
		pro_it = inMap.begin();

		for(; pro_it != inMap.end(); pro_it++){//inMap is the PPI network, size comparable to the size of sequence
			//compare pro_it->first with ai, if similar add (ai, pro_it->second) to nbMap
			pname=pro_it->first;
			//cout<<pname<<endl;
			seq="";

			if (seqmap.find(pname) != seqmap.end()){
				seq = seqmap[pname];
				//cout<<ai<<" "<<pname<<endl;
				//judge ai(a window) is similar to seq or not
				if( windows_similar(ai, seq)){
					add_to_NbMap(nbMap, ai, pro_it->second);
				}
			}//end if

		}//end for
		//cout<<nbMap[ai].size()<<endl;
	}//end outer for

	cout<<"Constructed Neighbor Map of A, size:"<<nbMap.size()<<endl;
	cout<<"----------------------------"<<endl;

	/*************************************************************************************/
    /****************************START TO PROCESS B***************************************/
	/*************************************************************************************/

	int abScore[bLen][aLen];  //store the score matrix of a_i, b_i
	for (int abi =0; abi<bLen; abi++)
	{
		for (int abj =0; abj<aLen; abj++)
		{
			abScore[abi][abj]=0;
		}
	}

	int max=0;//get the reliable number score of a_i, b_i
    ofstream maxoutput;
    maxoutput.open(SCORE_FILE.c_str());

	for(int i=0; i<=bLen-WINDOW_LEN; i++)
	{
		if (i==0){
			bi = b.substr(i, WINDOW_LEN);
		}
		else{
			bi = bi.erase(0, 1)+b[i+WINDOW_LEN-1];
		}

		//cout<<" bi:"<<bi<<endl;

		nbMap_it = nbMap.begin();
		int aindex=0;
		for(; nbMap_it != nbMap.end(); nbMap_it++)
		{
            vector<string> nbProtein= nbMap_it->second;
            //cout<<nbProtein.size()<<nbProtein[0]<<endl;
            for(unsigned int j=0; j<nbProtein.size(); j++)
            {
            	//cout<<" pname:"<<nbProtein[j]<<endl;
            	seq="";
    			if (seqmap.find( nbProtein[j] ) != seqmap.end())
    				seq = seqmap[nbProtein[j]];
    			//cout<<bi<<"\t"<<seq<<endl;
    			if( windows_similar(bi, seq))
    			{
    				//cout<<i<<" "<<aindex<<" "<<"similar"<<endl;
    				abScore[i][aindex]++;
    				//cout<<"score:"<<abScore[i][aindex]<<endl;
    			}

            }
            maxoutput<<abScore[i][aindex]<<" ";
            if (abScore[i][aindex]>max)
            {
            	max = abScore[i][aindex];
            }
			aindex ++;
			//cout<<"-------------"<<endl;

		}
		maxoutput<<endl;
	}

	if (max>=MATCH)
	{
		cout<<"Score:"<<max<<" A and B interact"<<endl;
	}
	else
	{
		cout<<"Score:"<<max<<" No interaction"<<endl;
	}

	maxoutput.close();
	return 1;
}

/*
 * a_i is k-mer str, seq is a protein sequence
 * right after check
 */
int windows_similar(string ai, string seq)
{
	int seqLen = seq.length();
	double score=-100000;
	double temp=0;

	for(int i=0; i<=seqLen-WINDOW_LEN; i++)
	{
		temp=0.0;
		//cout<<" ";
		for (unsigned int w=0 ; w<WINDOW_LEN; w++){
			//cout<<seq[i+w];
			temp += getPAMScore( ai[w], seq[i+w]);
		}

		//cout<<" "<<temp<<endl;
		if (temp>score)
			score = temp;

		if (score>=SPAM){
			//cout<<" add"<<endl;
			return 1;
		}
		//cout<<" e:"<<i<<" "<<score<<endl;
	}
	//cout<<endl;
	return 0;
}

/*
 * if find ai is similar to a node in graph, add its neighbors to nbMap
 * check if duplciate existed before adding
 */
int add_to_NbMap( hash_map<string, vector<string>, str_hash > &nbMap, string ai, vector<string> &prosecond)
{
	if( nbMap.find(ai)!=nbMap.end())
	{

		vector<string> nb = nbMap[ai];

		for( unsigned int i=0; i<prosecond.size(); i++)
		{

			if ( find(nb.begin(), nb.end(), prosecond[i]) == nb.end())
			{
				nbMap[ai].push_back(prosecond[i]);
				//cout<<"size:"<<nbMap[ai].size()<<endl;
			}
		}
	}
	else
	{
		nbMap[ai] = prosecond;
	}

	return 1;
}
/*
 * store seq infor to map
 * key is pname, value is the amino acid sequence
 */
int getSeqMap(string seqfile, hash_map<string, string, str_hash> & seqMap)
{
	ifstream input;
	input.open(seqfile.c_str());
	if (!input.is_open())
	{
		cout<<seqfile<<" open error"<<endl;
		return 0;
	}

	string line="";
	string pname="";
	string pseq="";
	while( getline(input, line))
	{
		if ( line.find(">")==0 )
		{
			pname = line.substr(1);
		}
		else
		{
			seqMap[pname]=line;
		}
	}
	cout<<"Sequence file Size: "<<seqMap.size()<<endl;
	input.close();
	return 1;
}

