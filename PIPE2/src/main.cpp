/*
 * main.cpp
 *
 *  Created on: Sep 30, 2014
 *      Author: yqian33
 */

#include "pipe.h"
#include "helper.h"

int main(void)
{
	outUsedTime(0);

	//readBlosum80(BLOSUM_FILE);
	readPam120(PAM_FILE);
	//cout<<"score(A,A):"<<getBlosumScore('A','A')<<endl;
	cout<<"score(*,A):"<<getPAMScore('*','X')<<endl;
	cout<<"Load score matrix done...\n----------------------------"<<endl;

	hash_map<string, vector<string> , str_hash>  inMap;
	interaction_Map(PPI_FILE,inMap);
	//cout<<check_inmap("O13329", inMap)<<endl;
	cout<<"----------------------------"<<endl;

	hash_map<string, string, str_hash> seqmap;
	getSeqMap(SEQ_FILE, seqmap);
	cout<<"----------------------------"<<endl;

	string a=seqmap["P32381"];
	string b=seqmap["Q05785"];
    //a="MSYTDNPPQTKRALSLDDLV";
    //b=a;
	pipe(a, b, inMap, seqmap);
	outUsedTime(1);

	return 0;
}


