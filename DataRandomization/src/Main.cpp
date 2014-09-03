/*
 * Main.cpp
 *
 *  Created on: Aug 28, 2014
 *      Author: yqian33
 */

#include "def.h"
void forAlignMCL()
{
		string dir="/home/yqian33/workspace/PPIDataProcess/AddInteraction/src/yeast-test/";
		string infile=dir+"fly-yeast/yeast.DIP.nif";
		int fly_total=24220;
		int yeast_total=22377;

		int number=5000;
		stringstream ss;
		ss<<number;

		int iter=10;
		int i=0;
		while(i<iter)
		{
			stringstream s1;
			s1<<i;
			string outfile=dir+"yeast_remove/"+"rm"+ss.str()+"_"+s1.str()+".nif";
			remove_randomly(infile, yeast_total, number, outfile);
			i++;
		}
}
void forMaWISh()
{
		string dir="/home/yqian33/APPI/software/MAWISH/data/ppi/";
		string infile=dir+"1/4932.ppi";
		int fly_total=24220;
		int yeast_total=22377;

		int number=1000;
		stringstream ss;
		ss<<number;

		int iter=1;
		int i=0;
		while(i<iter)
		{
			stringstream s1;
			s1<<i;
			string outfile=dir+"4932.ppi";
			remove_randomly(infile, yeast_total, number, outfile);
			i++;
		}
}

void add_to_alignmcl()
{
	string outfile="/home/yqian33/workspace/PPIDataProcess/AddInteraction/src/add_interaction/yeast.DIP.new.add5000.nif";
	string infile="/home/yqian33/workspace/PPIDataProcess/AddInteraction/src/add_interaction/origin/yeast.DIP.new.nif";

	//string yeastfile="/home/yqian33/workspace/PPIDataProcess/AddInteraction/src/yeast-test/yeast.pred1.4_sort_clean.txt";
	//string flyfile="/home/yqian33/workspace/PPIDataProcess/AddInteraction/src/fly.DIP-inverse/fly_pred1.5_sort_clean.txt";
	string yeastfile="/home/yqian33/APPI/software/AlignMCL-1.2/test/AlignMCLdataset/ppi/add/yeast.i2d.dip.nif";
	string flyfile="/home/yqian33/APPI/software/AlignMCL-1.2/test/AlignMCLdataset/ppi/add/fly.i2d.dip.nif";

	addInteraction(5000, infile, yeastfile, outfile);
	cout<<"done"<<endl;
}
int main(void)
{
	//forMaWISh();
	add_to_alignmcl();

}


